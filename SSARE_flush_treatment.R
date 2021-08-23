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
flush<-read.csv("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\SSARE_Flush eggs and nymphs_DATA.csv")
date<-flush$?..Date
time<-as.factor(flush$Timepoint)
grove<-flush$Grove
rep<-as.factor(flush$Rep)
plot<-flush$Plot
row<-flush$Row
tree<-flush$Tree
tree2<-flush$T
set<-flush$Set
tx<-flush$Treatment
nymph<-flush$ACP_nymph
eggs<-flush$ACP_egg
clm<-flush$prop_CLM

ordered.time<- factor(c("1","2","3","4", "5","6","7","8","9")) 


#######################################################################################################

#Comparing ACP nymphs per treatment

y.control.g.nymph<-tapply(nymph[which(tx=="control.g")], time[which(tx=="control.g")], mean, na.rm=TRUE)
#y.control.g.nymph<- y.cont[ordered.time]
sd.control.g.nymph<-tapply(nymph[which(tx=="control.g")], time[which(tx=="control.g")], sd, na.rm=TRUE)
#sd.control.g.nymph<-sd.cont[ordered.time]
n.control.g.nymph<-tapply(nymph[which(tx=="control.g")], time[which(tx=="control.g")], length)
#n.control.g.nymph<-n.cont[ordered.time]
se.control.g.nymph<-sd.control.g.nymph/sqrt(n.control.g.nymph)

y.control.nymph<-tapply(nymph[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.nymph<-tapply(nymph[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.nymph<-tapply(nymph[which(tx=="control")], time[which(tx=="control")], length)
se.control.nymph<-sd.control.nymph/sqrt(n.control.nymph)

y.flower.nymph<-tapply(nymph[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.nymph<-tapply(nymph[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.nymph<-tapply(nymph[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.nymph<-sd.flower.nymph/sqrt(n.flower.nymph)

y.vine.nymph<-tapply(nymph[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.nymph<-tapply(nymph[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.nymph<-tapply(nymph[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.nymph<-sd.vine.nymph/sqrt(n.vine.nymph)

y.bush.nymph<-tapply(nymph[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.nymph<-tapply(nymph[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.nymph<-tapply(nymph[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.nymph<-sd.bush.nymph/sqrt(n.bush.nymph)

y.all.nymph<-tapply(nymph[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.nymph<-tapply(nymph[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.nymph<-tapply(nymph[which(tx=="all")], time[which(tx=="all")], length)
se.all.nymph<-sd.all.nymph/sqrt(n.all.nymph)

#######################################################################################################

#Line graphs for average visitor nymphdance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_nymphs_tx_aug 17.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,150), type="n", xlab="Sampling time", ylab="Avg. nymph abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July","Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.g.nymph-se.control.g.nymph[1:length(unique(time))], 1:length(unique(time)), y.control.g.nymph+se.control.g.nymph[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.nymph-se.control.nymph[1:length(unique(time))], 1:length(unique(time)), y.control.nymph+se.control.nymph[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.nymph-se.flower.nymph[1:length(unique(time))], 1:length(unique(time)), y.flower.nymph+se.flower.nymph[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.nymph-se.vine.nymph[1:length(unique(time))], 1:length(unique(time)), y.vine.nymph+se.vine.nymph[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.nymph-se.bush.nymph[1:length(unique(time))], 1:length(unique(time)), y.bush.nymph+se.bush.nymph[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.nymph-se.all.nymph[1:length(unique(time))], 1:length(unique(time)), y.all.nymph+se.all.nymph[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control.g.nymph, type="b", col="black", bg="black",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.nymph, type="b", col="grey", bg="grey", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.nymph, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.nymph, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.nymph-0.025, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.nymph-0.05, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=150, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"), 
       col=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"))

dev.off()


#######################################################################################################

#Comparing eggs per treatment

y.control.g.eggs<-tapply(eggs[which(tx=="control.g")], time[which(tx=="control.g")], mean, na.rm=TRUE)
#y.control.g.eggs<- y.cont[ordered.time]
sd.control.g.eggs<-tapply(eggs[which(tx=="control.g")], time[which(tx=="control.g")], sd, na.rm=TRUE)
#sd.control.g.eggs<-sd.cont[ordered.time]
n.control.g.eggs<-tapply(eggs[which(tx=="control.g")], time[which(tx=="control.g")], length)
#n.control.g.eggs<-n.cont[ordered.time]
se.control.g.eggs<-sd.control.g.eggs/sqrt(n.control.g.eggs)

y.control.eggs<-tapply(eggs[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.eggs<-tapply(eggs[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.eggs<-tapply(eggs[which(tx=="control")], time[which(tx=="control")], length)
se.control.eggs<-sd.control.eggs/sqrt(n.control.eggs)

y.flower.eggs<-tapply(eggs[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.eggs<-tapply(eggs[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.eggs<-tapply(eggs[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.eggs<-sd.flower.eggs/sqrt(n.flower.eggs)

y.vine.eggs<-tapply(eggs[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.eggs<-tapply(eggs[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.eggs<-tapply(eggs[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.eggs<-sd.vine.eggs/sqrt(n.vine.eggs)

y.bush.eggs<-tapply(eggs[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.eggs<-tapply(eggs[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.eggs<-tapply(eggs[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.eggs<-sd.bush.eggs/sqrt(n.bush.eggs)

y.all.eggs<-tapply(eggs[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.eggs<-tapply(eggs[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.eggs<-tapply(eggs[which(tx=="all")], time[which(tx=="all")], length)
se.all.eggs<-sd.all.eggs/sqrt(n.all.eggs)

#######################################################################################################

#Line graphs for average eggs at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_eggs_tx_aug17.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,100), type="n", xlab="Sampling time", ylab="Avg. egg abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.g.eggs-se.control.g.eggs[1:length(unique(time))], 1:length(unique(time)), y.control.g.eggs+se.control.g.eggs[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.eggs-se.control.eggs[1:length(unique(time))], 1:length(unique(time)), y.control.eggs+se.control.eggs[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.eggs-se.flower.eggs[1:length(unique(time))], 1:length(unique(time)), y.flower.eggs+se.flower.eggs[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.eggs-se.vine.eggs[1:length(unique(time))], 1:length(unique(time)), y.vine.eggs+se.vine.eggs[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.eggs-se.bush.eggs[1:length(unique(time))], 1:length(unique(time)), y.bush.eggs+se.bush.eggs[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.eggs-se.all.eggs[1:length(unique(time))], 1:length(unique(time)), y.all.eggs+se.all.eggs[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control.g.eggs, type="b", col="black", bg="black",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.eggs, type="b", col="grey", bg="grey", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.eggs, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.eggs, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.eggs-0.025, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.eggs-0.05, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=100, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"), 
       col=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"))

dev.off()