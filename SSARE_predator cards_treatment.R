################################################################################################################
### Project: SSARE Windbreaks ##################################################################################
## Date: July 2021  ## Authors: Chuang, ... ## Year: 2021 ## Version: Meeting with Xavier  #####################
################################################################################################################
### Purpose: Analyze Flush Data by Grove by Timepoint ##########################################################
################################################################################################################

### Install packages, if needed 
if (!require(grDevices)) {  
  install.packages("grDevices", repos = "http://cran.us.r-project.org")
  require(grDevices)
}

if (!require(MASS)) {  
  install.packages("MASS", repos = "http://cran.us.r-project.org")
  require(MASS)
}

if (!require(car)) {  
  install.packages("car", repos = "http://cran.us.r-project.org")
  require(car)
}

if (!require(glmmTMB)) {  
  install.packages("glmmTMB", repos = "http://cran.us.r-project.org")
  require(glmmTMB)
}

if (!require(DHARMa)) {  
  install.packages("DHARMa", repos = "http://cran.us.r-project.org")
  require(DHARMa)
}

if (!require(chron)) {  
  install.packages("chron", repos = "http://cran.us.r-project.org")
  require(chron)
}

if (!require(vioplot)) {  
  install.packages("chron", repos = "http://cran.us.r-project.org")
  require(chron)
}

if (!require(emmeans)) {  
  install.packages("emmeans", repos = "http://cran.us.r-project.org")
  require(emmeans)
}

if (!require(multcomp)) {  
  install.packages("multcomp", repos = "http://cran.us.r-project.org")
  require(multcomp)
}

if (!require(wesanderson)) {  
  install.packages("wesanderson", repos = "http://cran.us.r-project.org")
  require(wesanderson)
}

if (!require(fitdistrplus)) {  
  install.packages("fitdistrplus", repos = "http://cran.us.r-project.org")
  require(fitdistrplus)
}

if (!require(scales)) {  
  install.packages("scales", repos = "http://cran.us.r-project.org")
  require(scales)
}

### Load libraries 
library(grDevices) #to save to jpeg
library(MASS)
library(car)
library(glmmTMB)
library(fitur)
library(DHARMa)
library(chron) #to make time objects
library(vioplot)
library(emmeans) #to run pairwise comparisons 
library(multcomp) #to set up pairwise comparisons and assign letters via function cld
library(wesanderson) #for pleasing color palette
library(fitdistrplus)
library(scales) #to use alpha for transparent colors

###Set up my working directory so the csv ends up in the right folder
setwd("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\")

######################
### Load dataset
pred<-read.csv("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\SSARE_Predator card DATA.csv")
date<-pred$?..Date
time<-as.factor(pred$Timepoint)
grove<-pred$Grove
rep<-as.factor(pred$Rep)
set<-pred$Set
plot<-pred$Plot
row<-pred$Row
tree<-pred$Tree
bag<-pred$Bagged
tx<-pred$Treatment
intact<-pred$Intact.prop

#######################################################################################################

#Comparing intact eggs per treatment

y.control.g.intact<-tapply(intact[which(tx=="control.g")], time[which(tx=="control.g")], mean, na.rm=TRUE)
#y.control.g.intact<- y.cont[ordered.time]
sd.control.g.intact<-tapply(intact[which(tx=="control.g")], time[which(tx=="control.g")], sd, na.rm=TRUE)
#sd.control.g.intact<-sd.cont[ordered.time]
n.control.g.intact<-tapply(intact[which(tx=="control.g")], time[which(tx=="control.g")], length)
#n.control.g.intact<-n.cont[ordered.time]
se.control.g.intact<-sd.control.g.intact/sqrt(n.control.g.intact)

y.control.intact<-tapply(intact[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.intact<-tapply(intact[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.intact<-tapply(intact[which(tx=="control")], time[which(tx=="control")], length)
se.control.intact<-sd.control.intact/sqrt(n.control.intact)

y.flower.intact<-tapply(intact[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.intact<-tapply(intact[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.intact<-tapply(intact[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.intact<-sd.flower.intact/sqrt(n.flower.intact)

y.vine.intact<-tapply(intact[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.intact<-tapply(intact[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.intact<-tapply(intact[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.intact<-sd.vine.intact/sqrt(n.vine.intact)

y.bush.intact<-tapply(intact[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.intact<-tapply(intact[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.intact<-tapply(intact[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.intact<-sd.bush.intact/sqrt(n.bush.intact)

y.all.intact<-tapply(intact[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.intact<-tapply(intact[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.intact<-tapply(intact[which(tx=="all")], time[which(tx=="all")], length)
se.all.intact<-sd.all.intact/sqrt(n.all.intact)

#######################################################################################################

#Line graphs for average visitor intactdance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\predator cards_intact_tx_jul21.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,2), type="n", xlab="Sampling time", ylab="Avg. proportion of intact eggs/sample", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.g.intact-se.control.g.intact[1:length(unique(time))], 1:length(unique(time)), y.control.g.intact+se.control.g.intact[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.intact-se.control.intact[1:length(unique(time))], 1:length(unique(time)), y.control.intact+se.control.intact[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.intact-se.flower.intact[1:length(unique(time))], 1:length(unique(time)), y.flower.intact+se.flower.intact[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.intact-se.vine.intact[1:length(unique(time))], 1:length(unique(time)), y.vine.intact+se.vine.intact[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.intact-se.bush.intact[1:length(unique(time))], 1:length(unique(time)), y.bush.intact+se.bush.intact[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.intact-se.all.intact[1:length(unique(time))], 1:length(unique(time)), y.all.intact+se.all.intact[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control.g.intact, type="b", col="black", bg="black",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.intact, type="b", col="grey", bg="grey", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.intact, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.intact, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.intact, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.intact, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=3, y=2, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"), 
       col=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"))

dev.off()


#######################################################################################################

