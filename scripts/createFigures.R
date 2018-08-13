#Load in the libraries we will need
library(ggplot2); library(rstan); library(metafor)
library(cowplot); library(grid); library(shinystan)

#Move into script directory
setwd(dirname(sys.frame(1)$ofile))

#Multicore stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Plotting options for ggplot2
theme_set(theme_gray())
ggOptions = theme(axis.text.x=element_text(size=10,family="Arial"),axis.text.y=element_text(size=10,family="Arial"),
                  axis.title.x=element_text(size=10,face="bold",family="Arial"),axis.title.y=element_text(size=10,angle=90,face="bold",family="Arial"),
                  legend.text=element_text(size=10,face="bold",family="Arial"),plot.title=element_text(size=12,face="bold",hjust=0.5,family="Arial"),
                  legend.title.align=0.5,legend.title=element_text(size=10,face="bold",family="Arial"))

#Load in  meta analysis data
metData = read.csv('../data/metaAvData.csv')

#Which measurment are we doing?
meas = "OCI"

#Make a generic variable for measurment variable
metData$Meas = get(meas,metData)
metData$Meas.Error = get(paste(meas,".Error",sep=''),metData)
metData$Meas.Error.Type = get(paste(meas,".Error.Type",sep=''),metData)
metData$Meas.N = get(paste(meas,".N",sep=''),metData)

#Remove nans
measData = metData[is.na(metData$Meas)==FALSE,]

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
  
  #Add study name and year
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

#Sort by year then alphabetically
refFactor = factor(uRef,levels=uRef[order(uYear,uRef,decreasing=c(TRUE,FALSE))],ordered=TRUE)

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
real<lower=0> theta;
} 

model {
muG ~ normal(6,2);
muS ~ normal(0,theta);
for (i in 1:N){
y[i] ~ normal(muG + muS[i], se[i]);
}
}

"

#Which parameters to save
params = c("muG","muS","theta")

#Make data list for stan
stanData = list(
  y = measFrame$meas,
  se = measFrame$se,
  N = length(measFrame$meas)
)

#Run stan
run = FALSE
if (run == TRUE){
  stanFit = stan(model_code=stanModel,data=stanData,pars=params,iter=20000,warmup=10000,chains=8,thin=5,control=list(adapt_delta=0.99))
  
  #Save posterior
  save("stanFit",file=paste('../data/',meas,'_posteriorSamples_',format(Sys.Date(),format='%m_%d_%y'),'.Rdata',sep=''))
}else{
  load(paste('../data/',meas,'_posteriorSamples_07_05_18.Rdata',sep=''))
}

#Get population parameters
muG = extract(stanFit,pars="muG")$muG
theta = extract(stanFit,pars="theta")$theta

#Compute within subjects variance
k = stanData$N
w = 1 / stanData$se^2
sigmaHat = sum(w)*(k-1) / (sum(w)^2 - sum(w^2))
bayesI = 100 * theta^2 / (theta^2+sigmaHat)

#Get quantiles
quantMuG = quantile(muG,probs=c(0.025,0.5,0.975))
quantTheta = quantile(theta,probs=c(0.025,0.5,0.975))
quantI = quantile(bayesI,probs=c(0.025,0.5,0.975))

#Show quantiles
print(paste(meas,' Population Mean:',sep='')); print(quantMuG)
print(paste(meas,' Theta:',sep='')); print(quantTheta)
print(paste(meas,' I^2:',sep='')); print(quantI)

#Get prediction interval
nSamples = dim(muG)
muPred = rep(0,nSamples)
for (i in 1:nSamples){
  muPred[i] = muG[i] + rnorm(1,0,theta[i])
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
                       study=factor(c(uRef,"","Overall"),levels=c(uRef,"","Overall")[order(c(uYear,2999,3000),c(uRef,'z','z'),decreasing=c(TRUE,FALSE))],ordered=TRUE))

#Make a combined plot
combPlot = ggplot(combFrame,aes(y=study,x=med))
combPoint = geom_point(color="#006BB6",fill="#006BB6",
                       size=combFrame$size,
                       shape=combFrame$shape)
combError = geom_segment(aes(x=low,xend=up,y=study,yend=study),size=1)
combLine = geom_vline(aes(xintercept=6.0),linetype="dashed",size=1)
combG = geom_hline(aes(yintercept=2.0),linetype="solid",size=1)
combX = scale_x_continuous(meas,limits=c(2,14.75),
                           breaks=seq(from=2,to=14,by=2),
                           labels=seq(from=2,to=14,by=2))
combY = scale_y_discrete("Study",drop=FALSE)
combTitle = ggtitle(paste(meas,'Meta-analysis'))
combFig = combPlot + combPoint + combLine + combError + combG
combFig = combFig + combTitle + combX + combY + ggOptions + theme(plot.margin = unit(c(1,9,1,1), "lines"))

nStudy = length(combFrame$med)
for (i in 1:nStudy)  {
  if (!is.na(combFrame$med[i])){
    combFig = combFig + annotation_custom(
      grob = textGrob(label = sprintf("%.2f [%.2f, %.2f]",
                                      combFrame$med[i],
                                      combFrame$low[i],combFrame$up[i]),
                      hjust = 0, gp = gpar(cex = 1.15,fontface=2)),
      ymin = which(levels(combFrame$study)==combFrame$study[i]),
      ymax = which(levels(combFrame$study)==combFrame$study[i]),
      xmin = 15.75,
      xmax = 15.75)
  } 
}  

#Override clipping
combTable = ggplot_gtable(ggplot_build(combFig))
combTable$layout$clip[combTable$layout$name == "panel"] = "off"

#Save plot
if (meas == "OGI"){
  outName = "Fig2.tif"
}else{
  outName = "Fig3.tif"
}
tiff(paste('../figures/',outName,sep=''),width=7,7,units='in',res=300)
grid.draw(combTable)
dev.off()

#Make dataframe for polygon
cFrame = data.frame(x=quantMuG,y=c(3.25,0,3.25))
pFrame = data.frame(x=quantPred,y=c(3.25,0,3.25))

#Make a quick funnel plot
fPlot = ggplot(measFrame,aes(x=meas,y=se))
fPoint = geom_point(size=4,color="#006BB6")
fY = scale_y_reverse("Standard Error",limits=c(3.25,0))
fX = scale_x_continuous(meas,limits=c(3,9))
fEst = geom_vline(xintercept=quantMuG[2],size=1,linetype="dashed")
fPolyC = geom_polygon(data=cFrame,aes(x=x,y=y),fill="gray65",color="black")
fPolyP = geom_polygon(data=pFrame,aes(x=x,y=y),fill="gray85",color="black")
fTitle = ggtitle(paste(meas,'Funnel Plot'))
fFigure = fPlot + fPolyP + fPolyC + fPoint + fX + fY + fEst + fTitle + ggOptions

#Save funnel plot
tiff(paste('../figures/',meas,'FunnelStudy.tiff',sep=''),width=6,height=6,units='in',res=300)
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
                       labels=c('A)','B)'),label_size=11,vjust=1,hjust=0)
  
  #Add title
  plotTitle = ggdraw() + draw_label("Funnel Plots", fontface='bold',size=12,fontfamily='Arial')
  plotComb = plot_grid(plotTitle, plotComb, ncol=1, rel_heights=c(0.1, 1))
  
  #Save final plot
  tiff('../figures/Fig4.tif',width=7,height=3.5,units='in',res=300)
  print(plotComb)
  dev.off()
  
}

#Run freq. meta analysis.
freqFit = rma(data = measFrame, yi = meas, sei = se, slab = ref)

#Run regression test
regTest = regtest(freqFit,ret.fit=TRUE)

#Show user regression test p-value
print(paste(meas,' Asymmetry p-value = ',round(regTest$pval,4),sep=''))

#Make a shiny stan object of fit
sStan = as.shinystan(stanFit)

#Get rhat
rHat = retrieve(sStan,"rhat")

#Show rHat that is futherest from zero
print(max(abs(1.0 - rHat)))

