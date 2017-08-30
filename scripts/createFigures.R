#Load in the libraries we will need
library(ggplot2); library(rstan); library(metafor); library(cowplot)

#Move into script directory
setwd(dirname(sys.frame(1)$ofile))

#Multicore stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Plotting options for ggplot2
theme_set(theme_gray())
ggOptions = theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14,face="bold"),
                  axis.title.y=element_text(size=14,angle=90,face="bold"),legend.text=element_text(size=14,face="bold"),
                  plot.title=element_text(size=16,face="bold",hjust=0.5),legend.title.align=0.5,legend.title=element_text(size=14,face="bold"))

#Load in  meta analysis data
metData = read.csv('../data/cleanData.csv')

#Which measurment are we doing?
meas = "OCI"

#Make a generic variable for measurment variable
metData$Meas = get(meas,metData)
metData$Meas.Error = get(paste(meas,".Error",sep=''),metData)
metData$Meas.Error.Type = get(paste(meas,".Error.Type",sep=''),metData)
metData$Meas.N = get(paste(meas,".N",sep=''),metData)

#Remove nans
measData = metData[is.na(metData$Meas)==FALSE&metData$Meas!="ASK",]

#Make numeric values
measData$Meas = as.numeric(as.character(measData$Meas))
measData$Meas.Error = as.numeric(as.character(measData$Meas.Error))
measData$Meas.N = as.numeric(as.character(measData$Meas.N))
measData$Year = as.numeric(measData$Year)

#Compute SE
measData$Meas.SE =  measData$Meas.Error
measData$Meas.SE[measData$Meas.Error.Type=="SD"] = measData$Meas.Error[measData$Meas.Error.Type=="SD"] / sqrt(measData$Meas.N[measData$Meas.Error.Type=="SD"])

#Get IDs for unique studies
uStudy = unique(measData$Study)

#Make empty places holders for study averages
uRef = NULL; uMeas = NULL; uSe = NULL; uYear = NULL

#Loop through each study
for (s in uStudy){
  
  #Get indicies  and number of experiments in current study
  sIdx = measData$Study==s
  nExp = sum(sIdx)
  
  #Add study name,year, and gender
  uRef = append(uRef,as.character(measData$Reference[sIdx][1]))
  uYear = append(uYear,measData$Year[sIdx][1])
  
  #If we only have one experiment, just directly use the data
  if (nExp == 1){
    uMeas = append(uMeas,measData$Meas[sIdx])
    uSe =append(uSe,measData$Meas.SE[sIdx])
  }
  #Otherwise, run a quick fixed effects meta analysis for that study
  else{
    fixMet = rma.uni(measData$Meas[sIdx],sei=measData$Meas.SE[sIdx],method="FE")
    uMeas = append(uMeas,fixMet$beta[1])
    uSe = append(uSe,fixMet$se)
  }
  
}

#Get suffix sorter
lChar = substr(uRef, nchar(uRef), nchar(uRef))
lChar[lChar!="c"&lChar!="b"&lChar!="a"] = 3
lChar[lChar=="c"] = 1; lChar[lChar=="b"] = 2; lChar[lChar=="a"] = 3

#Sort by year and then by suffix
refFactor = factor(uRef,levels=uRef[order(-uYear,lChar)],ordered=TRUE)

#Make a data frame of study averages
measFrame = data.frame(ref=refFactor,meas=uMeas,se=uSe,year=uYear)

#Define stan model for meta analysis
stanModel = "

  data {
    int<lower=0> N;
    vector[N] se;
    vector[N] y;
  }

  parameters {
    real muG;
    vector[N] muS;
    real<lower=0> sigma;
  } 

  model {
    muG ~ normal(6,2);
    muS ~ normal(0,sigma);
    for (i in 1:N){
      y[i] ~ normal(muG + muS[i], se[i]);
    }
}

"

#Which parameters to save
params = c("muG","muS","sigma")

#Make data list for stan
stanData = list(
  y = measFrame$meas,
  se = measFrame$se,
  N = length(measFrame$meas)
)

#Run stan
stanFit = stan(model_code=stanModel,data=stanData,pars=params,iter=20000,warmup=10000,chains=8,thin=5,control=list(adapt_delta=0.99))

#Save posterior
save("stanFit",file=paste('../data/',meas,'_posteriorSamples_',format(Sys.Date(),format='%m_%d_%y'),'.Rdata',sep=''))

#Get population parameters
muG = extract(stanFit,pars="muG")$muG
sigma = extract(stanFit,pars="sigma")$sigma

#Get populations quantiles
quantMuG = quantile(muG,probs=c(0.025,0.5,0.975))
quantSigma = quantile(sigma,probs=c(0.025,0.5,0.975))

#Show population quantiles
print(paste(meas,' Population Mean:',sep='')); print(quantMuG)
print(paste(meas,' Population Sigma:',sep='')); print(quantSigma)

#Get prediction interval
nSamples = dim(muG)
muPred = rep(0,nSamples)
for (i in 1:nSamples){
  muPred[i] = muG[i] + rnorm(1,0,sigma[i])
}

#Get quantiles for prediction interval
quantPred = quantile(muPred,probs=c(0.025,0.5,0.975))

#Show prediction intervals
print(paste(meas,' Prediction Intervals:',sep=''))
print(quantPred)

#Make a combined dataframe of studies and meta analysis
combFactors = append(append(as.character(measFrame$ref),""),"Overall")
combFrame = data.frame(med=c(measFrame$meas,NA,quantMuG[2]),
                       shape = c(rep(15,length(measFrame$meas)),NA,23),
                       size=c(rep(4,length(measFrame$meas)),NA,5),
                       type=c(rep("Literature",length(measFrame$meas)),NA,"Meta-analysis"),
                       low=c(measFrame$meas-measFrame$se*1.96,NA,quantMuG[1]),
                       up=c(measFrame$meas+measFrame$se*1.96,NA,quantMuG[3]),
                       study=factor(c(uRef,"","Overall"),levels=c(uRef,"","Overall")[order(-c(uYear,2999,3000),c(lChar,3,3))],ordered=TRUE))

#Make a combined plot
combPlot = ggplot(combFrame,aes(y=study,x=med))
combPoint = geom_point(color="#006BB6",fill="#006BB6",
                       size=combFrame$size,
                       shape=combFrame$shape)
combError = geom_segment(aes(x=low,xend=up,y=study,yend=study),size=1)
combLine = geom_vline(aes(xintercept=6.0),linetype="dashed",size=1)
combG = geom_hline(aes(yintercept=2.0),linetype="solid",size=1)
combText = annotate("text",x=9.75,y=combFrame$study,
                    label=ifelse(is.na(combFrame$med),"",sprintf("%.2f [%.2f, %.2f]",combFrame$med,combFrame$low,combFrame$up)),
                    parse=FALSE,size=4.5,fontface=2,na.rm=TRUE)
combX = scale_x_continuous(meas,limits=c(2.5,10.75))
combY = scale_y_discrete("Study",drop=FALSE)
combTitle = ggtitle(paste(meas,'Meta-analysis'))
combFig = combPlot + combError + combPoint + combLine + combG
combFig = combFig + combTitle + combX + combY + combText + ggOptions
print(combFig)

#Save plot
pdf(paste('../figures/',meas,'CombFigStudy.pdf',sep=''),width=9,height=9)
print(combFig)
dev.off()

#Make dataframe for polygon
cFrame = data.frame(x=quantMuG,y=c(1.25,0,1.25))
pFrame = data.frame(x=quantPred,y=c(1.25,0,1.25))

#Make a quick funnel plot
fPlot = ggplot(measFrame,aes(x=meas,y=se))
fPoint = geom_point(size=4,color="#006BB6")
fY = scale_y_reverse("Standard Error",limits=c(1.25,0))
fX = scale_x_continuous(meas,limits=c(3,8))
fEst = geom_vline(xintercept=quantMuG[2],size=1,linetype="dashed")
fPolyC = geom_polygon(data=cFrame,aes(x=x,y=y),fill="gray65",color="black")
fPolyP = geom_polygon(data=pFrame,aes(x=x,y=y),fill="gray85",color="black")
fTitle = ggtitle(paste(meas,'Funnel Plot'))
fFigure = fPlot + fPolyP + fPolyC + fPoint + fX + fY + fEst + fTitle + ggOptions
print(fFigure)

#Save funnel plot
pdf(paste('../figures/',meas,'FunnelStudy.pdf',sep=''),width=6,height=6)
print(fFigure)
dev.off()

#Save funnel so we can combine. Horrible hack
if (meas == "OGI"){
  ogiFunnel = fFigure 
} else {
  ociFunnel = fFigure
}

#If both plots exist make a combined funnel plot
if ( exists('ogiFunnel') & exists('ociFunnel')  ) {
  
  #Combine funnel plots
  plotComb = plot_grid(ogiFunnel+theme(plot.title=element_blank()),
                       ociFunnel+theme(plot.title=element_blank()),
                       labels=c('A)','B)'),label_size=16,vjust=1,hjust=0)
  
  #Add title
  plotTitle = ggdraw() + draw_label("Funnel Plots", fontface='bold',size=18)
  plotComb = plot_grid(plotTitle, plotComb, ncol=1, rel_heights=c(0.1, 1))
  print(plotComb)
  
  #Save final plot
  pdf('../figures/funnelPlots.pdf',width=10,height=5)
  print(plotComb)
  dev.off()
  
}

#Run freq. meta analysis to get I^2
freqFit = rma(data = measFrame, yi = meas, sei = se, slab = ref)

#Show user I^2
print(paste(meas," I^2 = ",round(freqFit$I2,2),'%',sep=''))

#Run regression test
regTest = regtest(freqFit,ret.fit=TRUE)

#Show user regression test p-value
print(paste(meas,' Asymmetry p-value = ',round(regTest$pval,4),sep=''))



