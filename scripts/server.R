#Load in R libraries
library(ggplot2); library(shiny); library(rstan); library(metafor)

#Little function for computing SE
errorCalc = function(eVal,eN,eType){
  eCalc = rep(0,length(eVal))
  for (i in 1:length(eVal)){
    if (is.na(eVal[i])){
      eCalc[i] = NA
    }
    else if (eType[i]=="SE"){
      eCalc[i] = eVal[i]
    }
    else{
      eCalc[i] = eVal[i] / sqrt(eN[i])
    }
  }
  return(eCalc)
}

#Plotting options for ggplot2
ggOptions = theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14,face="bold"),
                  axis.title.y=element_text(size=14,angle=90,face="bold"),legend.text=element_text(size=14,face="bold"),
                  plot.title=element_text(size=16,face="bold",hjust=0.5),legend.title.align=0.5,legend.title=element_text(size=14,face="bold"))

#Load in  meta analysis data
metData = read.csv('../data/cleanData.csv')

#Replace ASK with NAN
metData[metData=="ASK"] = NA

#Make sure numeric data is numeric
metData$OGI = as.numeric(as.character(metData$OGI)); metData$OCI = as.numeric(as.character(metData$OCI))
metData$OGI.Error = as.numeric(as.character(metData$OGI.Error)); metData$OCI.Error = as.numeric(as.character(metData$OCI.Error))
metData$OGI.N = as.numeric(as.character(metData$OGI.N)); metData$OCI.N = as.numeric(as.character(metData$OCI.N))

#Get IDs for unique studies
uStudy = unique(metData$Study)

#Make empty places holders for study averages
uRef = NULL; uOgi = NULL; uOgiSe = NULL; uOci = NULL; uOciSe = NULL; uYear = NULL

#Loop through unique studies
for (u in uStudy){
  
  #Get indicies  and number of experiments in current study
  sIdx = metData$Study==u
  nExp = sum(sIdx)
  
  #Add study name to list
  uRef = append(uRef,as.character(metData$Reference[sIdx][1]))
  uYear = append(uYear,metData$Year[sIdx][1])
  
  #If we only have one experiment, just directly use the data
  if (nExp == 1){
    uOgi = append(uOgi,metData$OGI[sIdx]); uOci = append(uOci,metData$OCI[sIdx])
    uOgiSe = append(uOgiSe,errorCalc(metData$OGI.Error[sIdx],metData$OGI.N[sIdx],metData$OGI.Error.Type[sIdx]))
    uOciSe = append(uOciSe,errorCalc(metData$OCI.Error[sIdx],metData$OCI.N[sIdx],metData$OCI.Error.Type[sIdx]))
  }
  #Otherwise, run a quick fixed effects meta analysis for that study
  else{
    if (sum(is.na(metData$OGI[sIdx])==FALSE)>0){
      ogiMet = rma.uni(na.omit(metData$OGI[sIdx]),
                       sei=errorCalc(na.omit(metData$OGI.Error[sIdx]),
                                     na.omit(metData$OGI.N[sIdx]),
                                     na.omit(metData$OGI.Error.Type[sIdx])),
                       method="FE")
      uOgi = append(uOgi,ogiMet$beta[1])
      uOgiSe = append(uOgiSe,ogiMet$se)
    }
    else{
      uOgi = append(uOgi,NA)
      uOgiSe = append(uOgiSe,NA)
    }
    if (sum(is.na(metData$OCI[sIdx])==FALSE)>0){
      ociMet = rma.uni(na.omit(metData$OCI[sIdx]),
                       sei=errorCalc(na.omit(metData$OCI.Error[sIdx]),
                                     na.omit(metData$OCI.N[sIdx]),
                                     na.omit(metData$OCI.Error.Type[sIdx])),
                       method="FE")
      uOci = append(uOci,ociMet$beta[1])
      uOciSe = append(uOciSe,ociMet$se)
    }
    else{
      uOci = append(uOci,NA)
      uOciSe = append(uOciSe,NA)
    }
  }
}

#Make ogi and oci data lists
ogiList = list(y = as.numeric(na.omit(uOgi)),se = as.numeric(na.omit(uOgiSe)),
               N = length(na.omit(uOgi)),ref = uRef[is.na(uOgi)==FALSE],
               year = uYear[is.na(uOgi)==FALSE])
ociList = list(y = as.numeric(na.omit(uOci)),se = as.numeric(na.omit(uOciSe)),
               N = length(na.omit(uOci)),ref = uRef[is.na(uOci)==FALSE],
               year = uYear[is.na(uOci)==FALSE])

#Sort lists
ogiOrder = order(ogiList$year); ociOrder = order(ociList$year);
ogiList$y = ogiList$y[ogiOrder]; ociList$y = ociList$y[ociOrder];
ogiList$se = ogiList$se[ogiOrder]; ociList$se = ociList$se[ociOrder];
ogiList$ref = ogiList$ref[ogiOrder]; ociList$ref = ociList$ref[ociOrder];
ogiList$year = ogiList$year[ogiOrder]; ociList$year =ociList$year[ociOrder];

#Remove nans
ogiList = na.omit(ogiList)
ociList = na.omit(ociList)

#Define Stan model
stanModel = "

  data {
    int<lower=0> N;
    vector[N] se;
    vector[N] y;
  }

  parameters {
    real<lower=0,upper=20> muG;
    vector[N] muS;
    real<lower=0> sigma;
  } 

  model {
    muS ~ normal(0,sigma);
    for (i in 1:N){
      y[i] ~ normal(muG + muS[i], se[i]);
    }
  }

"

#Which parameters to save
params = c("muG","muS","sigma")

#Server function for Shiny app
function(input,output,session) {
  
  #Set initial reactive values
  rValues = reactiveValues(dataList=ogiList,useStudy=rep(FALSE,ogiList$N),label="OGI",
                           muG=NULL,sigma=NULL,plotFrame=NULL)
  
  #Reactive function to set appropriate data list
  makeFrame = reactive({
    
    if (input$meas == 1){
      rValues$dataList = ogiList
      rValues$useStudy = rep(FALSE,ogiList$N)
      rValues$label = "OGI"
    }
    else{
      rValues$dataList = ociList
      rValues$useStudy = rep(FALSE,ociList$N)
      rValues$label = "OCI"
    }
    
  })
  
  #Reactive funciton to run stan
  runStan = reactive({
    
    #Update dataframe
    makeFrame()
    
    #Make stan data frame
    stanData = list(y=rValues$dataList$y[rValues$useStudy==FALSE],
                    se=rValues$dataList$se[rValues$useStudy==FALSE],
                    N=sum(rValues$useStudy==FALSE))
    
    #Run stan
    if (file.exists('../data/stanModel.Rdata')){
      load('../data/stanModel.Rdata')
    } else{
      stanFit = stan(model_code=stanModel,data=stanData,save_dso=TRUE,pars=params,iter=0,warmup=0,chains=0,control=list(adapt_delta=0.99))
      save('stanFit',file='../data/stanModel.Rdata')
    }
    stanFit = stan(model_code=stanModel,data=stanData,save_dso=TRUE,pars=params,iter=2500,warmup=1250,chains=2,control=list(adapt_delta=0.99))
    
    #Get group estimates
    rValues$muG = extract(stanFit,pars="muG")$muG
    rValues$sigma = extract(stanFit,pars="sigma")$sigma
    
  })
  
  #Make meta-analysis plot
  output$metaPlot = renderPlot({
    
    #Update stan if necessary
    runStan()
    
    #Get group quantitlies
    lower = (1 - input$intLength/100) / 2.0
    upper = 1 - lower
    gQuant = quantile(rValues$muG,probs=c(lower,0.5,upper))
    seMul = qnorm(upper)
    
    #Make a combined dataframe with group result
    plotList = rValues$dataList
    
    #Make factors for plot
    plotFactors = append(append(as.character(plotList$ref),""),"Overall")
    plotFactors = reorder(factor(plotFactors),-c(plotList$year,2999,3000))
    
    #Sort factor again by a/b
    lChar = substr(plotList$ref, nchar(plotList$ref), nchar(plotList$ref))
    lChar[lChar!="b"] = 2; lChar[lChar=="b"] = 1
    plotFactors = reorder(factor(plotFactors),as.numeric(c(lChar,2,2)))
    
    #Make colors for plot
    plotColors = rep("#006BB6",plotList$N+2)
    plotColors[rValues$useStudy] = "white"
    plotFrame = data.frame(med=c(plotList$y,NA,gQuant[2]),
                           shape = c(rep(15,plotList$N),NA,23),
                           type=c(rep("Literature",plotList$N),NA,"Meta-analysis"),
                           low=c(plotList$y-plotList$se*seMul,NA,gQuant[1]),
                           up=c(plotList$y+plotList$se*seMul,NA,gQuant[3]),
                           study=reorder(factor(plotFactors),-c(plotList$year,2999,3000)),
                           color=plotColors,
                           size=c(rep(4,plotList$N),NA,5))
    
    #Save dataframe
    rValues$plotFrame = plotFrame
   
    #Make a plot with combined data
    combPlot = ggplot(plotFrame,aes(y=study,x=med))
    combPoint = geom_point(color=plotFrame$color,fill="#006BB6",
                           size=plotFrame$size,
                           shape=plotFrame$shape)
    combError = geom_segment(aes(x=low,xend=up,y=study,yend=study),size=1)
    combLine = geom_vline(aes(xintercept=6.0),linetype="dashed",size=1)
    combG = geom_hline(aes(yintercept=2.0),linetype="solid",size=1)
    combText = annotate("text",x=9.75,y=plotFrame$study,
                        label=ifelse(is.na(plotFrame$med),"",sprintf("%.2f [%.2f, %.2f]",plotFrame$med,plotFrame$low,plotFrame$up)),
                        parse=FALSE,size=4.5,fontface=2,na.rm=TRUE)
    combX = scale_x_continuous(rValues$label,limits=c(2.5,10.75))
    combY = scale_y_discrete("Study",drop=FALSE)
    combTitle = ggtitle(paste(rValues$label,"Meta-Analysis"))
    combFig = combPlot + combError + combPoint + combLine + combG
    combFig = combFig + combTitle + combX + combY + combText + ggOptions
    
    #Show figure
    makeOut()
    combFig
    
  })
  
  # Toggle points that are clicked
  observeEvent(input$metaClick, {
    
    #Get coordinate for study
    sStudy = (round(input$metaClick$y)-2)*-1 + rValues$dataList$N+1
    
    if (sStudy <= rValues$dataList$N+1){
      rValues$useStudy[sStudy] = !rValues$useStudy[sStudy]
    }
    
  })
  
  #Reactive function for cleaning up output for download
  makeOut = reactive({
    
    #Make output dataframe
    outFrame = rValues$plotFrame
    
    #Rename columns
    colnames(outFrame) = c(rValues$label,"Shape","Source",paste("CI",input$intLength,"Low",sep="."),
                           paste("CI",input$intLength,"High",sep="."),"Study","Include","Size")
    
    #Get rid of bad columns
    outFrame = na.omit(outFrame[,-c(2,8)])
    
    #Round to two digits
    outFrame[,1] = round(outFrame[,1],2)
    outFrame[,3] = round(outFrame[,3],2)
    outFrame[,4] = round(outFrame[,4],2)
    
    #Rename include variables
    outFrame$Include = as.character(outFrame$Include)
    outFrame$Include[outFrame$Include=="white"] = "No"
    outFrame$Include[outFrame$Include!="No"] = "Yes"
    
    #Sort columns
    outFrame = outFrame[,c(5,6,2,1,3,4)]
    
    #Return!
    return(outFrame)
    
  })
  
  #Downloader function
  output$download <- downloadHandler(
    
    #Write out data
    filename = function() { paste(rValues$label, '.csv', sep='') },
    content = function(file) {
      write.csv(makeOut(), file,row.names=FALSE)
    }
    
  )
  
  #Reset studies that are used
  observeEvent(input$reset, {
    rValues$useStudy = rep(FALSE,rValues$dataList$N)
  })
  
}

  