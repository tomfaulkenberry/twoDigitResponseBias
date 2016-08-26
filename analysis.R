########################################################################################################################
## R script for Faulkenberry, Cruise, & Shaki, "Reversing the manual digit bias in two-digit number comparison"
##
## note: there are helper functions at the bottom of the script that need to be excuted before plotting
########################################################################################################################

library(ggplot2)
rawData<-read.table("data.csv",sep=",",header=TRUE)

# clean up data
dataStep3<-subset(rawData,subset=error!=1) # remove errors

meanRT<-mean(dataStep3$RT)
sdRT<-sd(dataStep3$RT)
data<-subset(dataStep3,subset=RT<meanRT+3*sdRT & RT>200) # remove 3 SD outliers
attach(data)

meanInit<-mean(dataStep3$init.time)
sdInit<-sd(dataStep3$init.time)

meanMT<-mean(dataStep3$RT-dataStep3$init.time)
sdMT<-sd(dataStep3$RT-dataStep3$init.time)

########################################################################
## holistic processing signatures
########################################################################

# Time measures

# MT
agg=aggregate(RT-init.time~subject+distance+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
names(agg)[5]<-c("MT")
MT.aov=aov(MT~as.factor(distance)*side*mapping+Error(as.factor(subject)/(as.factor(distance)*side*mapping)),data=agg)
summary(MT.aov)
print(model.tables(MT.aov,"means"),digits=6)

# plot line graph (be sure to execute code at bottom to define "summarySEwithin" function)  Note: error bars defined in Morey 2008
summary=summarySEwithin(agg,measurevar="MT",withinvars=c("distance","mapping","side"),idvar="subject")
ggplot(summary,aes(x=distance,y=MT,shape=mapping))+geom_line(aes(group=mapping,linetype=mapping))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=MT-ci,ymax=MT+ci))+labs(x="Distance",y="Mean MT (ms)")+facet_grid(~side)+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))

# INIT
agg=aggregate(init.time~subject+distance+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
init.aov=aov(init.time~as.factor(distance)*side*mapping+Error(as.factor(subject)/(as.factor(distance)*side*mapping)),data=agg)
summary(init.aov)
print(model.tables(init.aov,"means"),digits=6)

#########################################################################################################
# plot hand trajectories for DISTANCE (close versus far)
# SE width at each timestep computed as standard error of mean x-coordinates over a sample of 64 participants
#########################################################################################################


dataLeftCloseCongruent<-subset(data,side=="left" & distance=="close" & mapping=="congruent")
dataRightCloseCongruent<-subset(data,side=="right" & distance=="close" & mapping=="congruent")
dataLeftFarCongruent<-subset(data,side=="left" & distance=="far" & mapping=="congruent")
dataRightFarCongruent<-subset(data,side=="right" & distance=="far" & mapping=="congruent")

dataLeftCloseIncongruent<-subset(data,side=="left" & distance=="close" & mapping=="incongruent")
dataRightCloseIncongruent<-subset(data,side=="right" & distance=="close" & mapping=="incongruent")
dataLeftFarIncongruent<-subset(data,side=="left" & distance=="far" & mapping=="incongruent")
dataRightFarIncongruent<-subset(data,side=="right" & distance=="far" & mapping=="incongruent")

# SE measures for each subset

SEmatrixLeftCloseCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftCloseIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftFarCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftFarIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightCloseCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightCloseIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightFarCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightFarIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)

for (i in 1:37){
  leftCloseCongruent<-subset(dataLeftCloseCongruent,subject==i)
  leftCloseIncongruent<-subset(dataLeftCloseIncongruent,subject==i)
  leftFarCongruent<-subset(dataLeftFarCongruent,subject==i)
  leftFarIncongruent<-subset(dataLeftFarIncongruent,subject==i)
  rightCloseCongruent<-subset(dataRightCloseCongruent,subject==i)
  rightCloseIncongruent<-subset(dataRightCloseIncongruent,subject==i)
  rightFarCongruent<-subset(dataRightFarCongruent,subject==i)
  rightFarIncongruent<-subset(dataRightFarIncongruent,subject==i)
  
  for (j in 1:101){
    SEmatrixLeftCloseCongruent[i,j]<-mean(leftCloseCongruent[,j+27])
    SEmatrixLeftCloseIncongruent[i,j]<-mean(leftCloseIncongruent[,j+27])
    SEmatrixLeftFarCongruent[i,j]<-mean(leftFarCongruent[,j+27])
    SEmatrixLeftFarIncongruent[i,j]<-mean(leftFarIncongruent[,j+27])
    SEmatrixRightCloseCongruent[i,j]<-mean(rightCloseCongruent[,j+27])
    SEmatrixRightCloseIncongruent[i,j]<-mean(rightCloseIncongruent[,j+27])
    SEmatrixRightFarCongruent[i,j]<-mean(rightFarCongruent[,j+27])
    SEmatrixRightFarIncongruent[i,j]<-mean(rightFarIncongruent[,j+27])
  }
}



xCoords=rep(0,808)
yCoords=rep(0,808)
decision=rep(0,808)
distance=rep(0,808)
mapping=rep(0,808)
SE=rep(0,808)

for (i in 1:101){
  xCoords[i]=mean(dataLeftCloseCongruent[,i+27])
  yCoords[i]=mean(dataLeftCloseCongruent[,i+128])
  decision[i]="smaller"
  distance[i]="close"
  mapping[i]="congruent"
  SE[i]=sd(SEmatrixLeftCloseCongruent[,i])/sqrt(37)

  xCoords[i+101]=mean(dataRightCloseCongruent[,i+27])
  yCoords[i+101]=mean(dataRightCloseCongruent[,i+128])
  decision[i+101]="larger"
  distance[i+101]="close"
  mapping[i+101]="congruent"
  SE[i+101]=sd(SEmatrixRightCloseCongruent[,i])/sqrt(37)
  
  xCoords[i+202]=mean(dataLeftFarCongruent[,i+27])
  yCoords[i+202]=mean(dataLeftFarCongruent[,i+128])
  decision[i+202]="smaller"
  distance[i+202]="far"
  mapping[i+202]="congruent"
  SE[i+202]=sd(SEmatrixLeftFarCongruent[,i])/sqrt(37)
  
  xCoords[i+303]=mean(dataRightFarCongruent[,i+27])
  yCoords[i+303]=mean(dataRightFarCongruent[,i+128])
  decision[i+303]="larger"
  distance[i+303]="far"
  mapping[i+303]="congruent"
  SE[i+303]=sd(SEmatrixRightFarCongruent[,i])/sqrt(37)
  
  
  xCoords[i+404]=mean(dataLeftCloseIncongruent[,i+27])
  yCoords[i+404]=mean(dataLeftCloseIncongruent[,i+128])
  decision[i+404]="larger"
  distance[i+404]="close"
  mapping[i+404]="incongruent"
  SE[i+404]=sd(SEmatrixLeftCloseIncongruent[,i])/sqrt(37)
  
  
  xCoords[i+505]=mean(dataRightCloseIncongruent[,i+27])
  yCoords[i+505]=mean(dataRightCloseIncongruent[,i+128])
  decision[i+505]="smaller"
  distance[i+505]="close"
  mapping[i+505]="incongruent"
  SE[i+505]=sd(SEmatrixRightCloseIncongruent[,i])/sqrt(37)
  
  xCoords[i+606]=mean(dataLeftFarIncongruent[,i+27])
  yCoords[i+606]=mean(dataLeftFarIncongruent[,i+128])
  decision[i+606]="larger"
  distance[i+606]="far"
  mapping[i+606]="incongruent"
  SE[i+606]=sd(SEmatrixLeftFarIncongruent[,i])/sqrt(37)
  
  
  xCoords[i+707]=mean(dataRightFarIncongruent[,i+27])
  yCoords[i+707]=mean(dataRightFarIncongruent[,i+128])
  decision[i+707]="smaller"
  distance[i+707]="far"
  mapping[i+707]="incongruent"
  SE[i+707]=sd(SEmatrixRightFarIncongruent[,i])/sqrt(37)
  
  
  
}

library("ggplot2")
trajectoryData=data.frame(xCoords,yCoords,decision,distance,mapping,SE)
plot=ggplot(trajectoryData,aes(x=yCoords,y=xCoords,group=mapping))+xlim(0,1.5)+ylim(-1,1)
paths=geom_path(aes(linetype=mapping),size=0.8)
ribbon=geom_ribbon(aes(x=yCoords,y=xCoords,ymin=xCoords-SE,ymax=xCoords+SE,linetype=mapping),alpha=0.2)

labels=labs(x="y-coordinates",y="x-coordinates")
faceting=facet_grid(decision~distance)
stripFormat=theme(strip.text=element_text(face="bold",size=rel(1.5)))
legendFormat=theme(legend.title=element_text(face="bold",size=rel(1.5)),legend.text=element_text(size=rel(1.5)))
axesFormat=theme(axis.title=element_text(size=rel(1.4)))
classic=theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))


basePlot=plot+paths+ribbon+labels+faceting+stripFormat+legendFormat+axesFormat+coord_flip()+classic
basePlot+labs(colour="Condition")+theme(legend.position=c(0.5,0.6))+theme(legend.background=element_rect(fill="white",colour="black"))

###############
## AUC analysis
###############

agg=aggregate(AUC~subject+distance+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
AUC.aov=aov(AUC~as.factor(distance)*side*mapping+Error(as.factor(subject)/(as.factor(distance)*side*mapping)),data=agg)
summary(AUC.aov)
print(model.tables(AUC.aov,"means"),digits=6)

# plot line graph (be sure to execute code at bottom to define "summarySEwithin" function)  Note: error bars defined in Morey 2008
summary=summarySEwithin(agg,measurevar="AUC",withinvars=c("distance","mapping","side"),idvar="subject")
ggplot(summary,aes(x=distance,y=AUC,shape=mapping))+geom_line(aes(group=mapping,linetype=mapping))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=AUC-ci,ymax=AUC+ci))+labs(x="Distance",y="Mean AUC")+facet_grid(~side)+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))




####################################################################################
##  Decomposed processing signatures
####################################################################################

## Time measures
# MT
agg=aggregate(RT-init.time~subject+compatibility+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
names(agg)[5]<-c("MT")
MT.aov=aov(MT~as.factor(compatibility)*side*mapping+Error(as.factor(subject)/(as.factor(compatibility)*side*mapping)),data=agg)
summary(MT.aov)
print(model.tables(MT.aov,"means"),digits=6)

# plot line graph (be sure to execute code at bottom to define "summarySEwithin" function)  Note: error bars defined in Morey 2008
summary=summarySEwithin(agg,measurevar="MT",withinvars=c("compatibility","mapping","side"),idvar="subject")
ggplot(summary,aes(x=compatibility,y=MT,shape=mapping))+geom_line(aes(group=mapping,linetype=mapping))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=MT-ci,ymax=MT+ci))+labs(x="Unit-decade compatibility",y="Mean MT (ms)")+facet_grid(~side)+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))
`
# INIT
agg=aggregate(init.time~subject+compatibility+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
init.aov=aov(init.time~as.factor(compatibility)*side*mapping+Error(as.factor(subject)/(as.factor(compatibility)*side*mapping)),data=agg)
summary(init.aov)
print(model.tables(init.aov,"means"),digits=6)


### Trajectory analyses

# plot hand trajectories for COMPATIBILITY

dataLeftCompatibleCongruent<-subset(data,side=="left" & compatibility=="compatible" & mapping=="congruent")
dataRightCompatibleCongruent<-subset(data,side=="right" & compatibility=="compatible" & mapping=="congruent")
dataLeftIncompatibleCongruent<-subset(data,side=="left" & compatibility=="incompatible" & mapping=="congruent")
dataRightIncompatibleCongruent<-subset(data,side=="right" & compatibility=="incompatible" & mapping=="congruent")

dataLeftCompatibleIncongruent<-subset(data,side=="left" & compatibility=="compatible" & mapping=="incongruent")
dataRightCompatibleIncongruent<-subset(data,side=="right" & compatibility=="compatible" & mapping=="incongruent")
dataLeftIncompatibleIncongruent<-subset(data,side=="left" & compatibility=="incompatible" & mapping=="incongruent")
dataRightIncompatibleIncongruent<-subset(data,side=="right" & compatibility=="incompatible" & mapping=="incongruent")

# SE measures for each subset

SEmatrixLeftCompatibleCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftCompatibleIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftIncompatibleCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixLeftIncompatibleIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightCompatibleCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightCompatibleIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightIncompatibleCongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)
SEmatrixRightIncompatibleIncongruent<-matrix(rep(0,37*101),nrow=37,ncol=101,byrow=TRUE)

for (i in 1:37){
  leftCompatibleCongruent<-subset(dataLeftCompatibleCongruent,subject==i)
  leftCompatibleIncongruent<-subset(dataLeftCompatibleIncongruent,subject==i)
  leftIncompatibleCongruent<-subset(dataLeftIncompatibleCongruent,subject==i)
  leftIncompatibleIncongruent<-subset(dataLeftIncompatibleIncongruent,subject==i)
  rightCompatibleCongruent<-subset(dataRightCompatibleCongruent,subject==i)
  rightCompatibleIncongruent<-subset(dataRightCompatibleIncongruent,subject==i)
  rightIncompatibleCongruent<-subset(dataRightIncompatibleCongruent,subject==i)
  rightIncompatibleIncongruent<-subset(dataRightIncompatibleIncongruent,subject==i)
  
  for (j in 1:101){
    SEmatrixLeftCompatibleCongruent[i,j]<-mean(leftCompatibleCongruent[,j+27])
    SEmatrixLeftCompatibleIncongruent[i,j]<-mean(leftCompatibleIncongruent[,j+27])
    SEmatrixLeftIncompatibleCongruent[i,j]<-mean(leftIncompatibleCongruent[,j+27])
    SEmatrixLeftIncompatibleIncongruent[i,j]<-mean(leftIncompatibleIncongruent[,j+27])
    SEmatrixRightCompatibleCongruent[i,j]<-mean(rightCompatibleCongruent[,j+27])
    SEmatrixRightCompatibleIncongruent[i,j]<-mean(rightCompatibleIncongruent[,j+27])
    SEmatrixRightIncompatibleCongruent[i,j]<-mean(rightIncompatibleCongruent[,j+27])
    SEmatrixRightIncompatibleIncongruent[i,j]<-mean(rightIncompatibleIncongruent[,j+27])
  }
}

xCoords=rep(0,808)
yCoords=rep(0,808)
side=rep(0,808)
compatibility=rep(0,808)
mapping=rep(0,808)
SE=rep(0,808)

for (i in 1:101){
  xCoords[i]=mean(dataLeftCompatibleCongruent[,i+27])
  yCoords[i]=mean(dataLeftCompatibleCongruent[,i+128])
  side[i]="left"
  compatibility[i]="compatible"
  mapping[i]="congruent"
  SE[i]=sd(SEmatrixLeftCompatibleCongruent[,i])/sqrt(37)
  
  
  xCoords[i+101]=mean(dataRightCompatibleCongruent[,i+27])
  yCoords[i+101]=mean(dataRightCompatibleCongruent[,i+128])
  side[i+101]="right"
  compatibility[i+101]="compatible"
  mapping[i+101]="congruent"
  SE[i+101]=sd(SEmatrixRightCompatibleCongruent[,i])/sqrt(37)
  
  
  xCoords[i+202]=mean(dataLeftIncompatibleCongruent[,i+27])
  yCoords[i+202]=mean(dataLeftIncompatibleCongruent[,i+128])
  side[i+202]="left"
  compatibility[i+202]="incompatible"
  mapping[i+202]="congruent"
  SE[i+202]=sd(SEmatrixLeftIncompatibleCongruent[,i])/sqrt(37)
  
  
  xCoords[i+303]=mean(dataRightIncompatibleCongruent[,i+27])
  yCoords[i+303]=mean(dataRightIncompatibleCongruent[,i+128])
  side[i+303]="right"
  compatibility[i+303]="incompatible"
  mapping[i+303]="congruent"
  SE[i+303]=sd(SEmatrixRightIncompatibleCongruent[,i])/sqrt(37)
  
  
  xCoords[i+404]=mean(dataLeftCompatibleIncongruent[,i+27])
  yCoords[i+404]=mean(dataLeftCompatibleIncongruent[,i+128])
  side[i+404]="left"
  compatibility[i+404]="compatible"
  mapping[i+404]="incongruent"
  SE[i+404]=sd(SEmatrixLeftCompatibleIncongruent[,i])/sqrt(37)
  
  
  xCoords[i+505]=mean(dataRightCompatibleIncongruent[,i+27])
  yCoords[i+505]=mean(dataRightCompatibleIncongruent[,i+128])
  side[i+505]="right"
  compatibility[i+505]="compatible"
  mapping[i+505]="incongruent"
  SE[i+505]=sd(SEmatrixRightCompatibleIncongruent[,i])/sqrt(37)
  
  
  xCoords[i+606]=mean(dataLeftIncompatibleIncongruent[,i+27])
  yCoords[i+606]=mean(dataLeftIncompatibleIncongruent[,i+128])
  side[i+606]="left"
  compatibility[i+606]="incompatible"
  mapping[i+606]="incongruent"
  SE[i+606]=sd(SEmatrixLeftIncompatibleIncongruent[,i])/sqrt(37)
  
  
  xCoords[i+707]=mean(dataRightIncompatibleIncongruent[,i+27])
  yCoords[i+707]=mean(dataRightIncompatibleIncongruent[,i+128])
  side[i+707]="right"
  compatibility[i+707]="incompatible"
  mapping[i+707]="incongruent"
  SE[i+707]=sd(SEmatrixRightIncompatibleIncongruent[,i])/sqrt(37)
  
  
  
}

library(ggplot2)
trajectoryData=data.frame(xCoords,yCoords,side,compatibility,mapping,SE)
plot=ggplot(trajectoryData,aes(x=yCoords,y=xCoords,group=compatibility))+xlim(0,1.5)+ylim(-1,1)
paths=geom_path(aes(linetype=compatibility),size=0.8)
ribbon=geom_ribbon(aes(x=yCoords,y=xCoords,ymin=xCoords-SE,ymax=xCoords+SE,linetype=compatibility),alpha=0.2)
labels=labs(x="y-coordinates",y="x-coordinates")
faceting=facet_grid(mapping~side)
stripFormat=theme(strip.text=element_text(face="bold",size=rel(1.5)))
legendFormat=theme(legend.title=element_text(face="bold",size=rel(1.5)),legend.text=element_text(size=rel(1.5)))
axesFormat=theme(axis.title=element_text(size=rel(1.4)))
classic=theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))


basePlot=plot+paths+ribbon+labels+faceting+stripFormat+legendFormat+axesFormat+coord_flip()+classic
basePlot+labs(colour="Condition")+theme(legend.position=c(0.5,0.6))+theme(legend.background=element_rect(fill="white",colour="black"))

### AUC analysis

agg=aggregate(AUC~subject+compatibility+side+mapping,data=data,FUN="mean") # RT performance data aggregated by subject
AUC.aov=aov(AUC~as.factor(compatibility)*side*mapping+Error(as.factor(subject)/(as.factor(compatibility)*side*mapping)),data=agg)
summary(AUC.aov)
print(model.tables(AUC.aov,"means"),digits=6)

# plot line graph (be sure to execute code at bottom to define "summarySEwithin" function)  Note: error bars defined in Morey 2008
summary=summarySEwithin(agg,measurevar="AUC",withinvars=c("compatibility","mapping","side"),idvar="subject")
ggplot(summary,aes(x=compatibility,y=AUC,shape=mapping))+geom_line(aes(group=mapping,linetype=mapping))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=AUC-ci,ymax=AUC+ci))+labs(x="Unit-decade compatibility",y="Mean AUC")+facet_grid(~side)+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))








# helper function for error bars above

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  require(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
