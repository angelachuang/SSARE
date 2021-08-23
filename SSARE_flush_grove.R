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


# Calling rm() Function to remove all objects
rm(list = ls())

######################
### Load dataset
flush<-read.csv("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\SSARE_Flush eggs and nymphs_DATA.csv")
date<-flush$?..Date
time<-as.factor(flush$Timepoint)
grove<-factor(flush$Grove, levels = c("control","enhanced"))
rep<-as.factor(flush$Rep)
plot<-flush$Plot
row<-flush$Row
tree<-flush$Tree
tree2<-flush$T
set<-flush$Set
tx<-flush$Treatment
nymph<-flush$ACP_nymph
eggs<-flush$ACP_egg
clm<-as.numeric(flush$prop_CLM)

ordered.time<- factor(c("1","2","3","4", "5","6","7","8","9")) 

#Comparing ACP nymphs in the enhanced vs. control grove 

y.control.n<-tapply(nymph[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.n<- y.cont[ordered.time]
sd.control.n<-tapply(nymph[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.n<-sd.cont[ordered.time]
n.control.n<-tapply(nymph[which(grove== "control")], time[which(grove== "control")], length)
#n.control.n<-n.cont[ordered.time]
se.control.n<-sd.control.n/sqrt(n.control.n)

y.enhanced.n<-tapply(nymph[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.n<-tapply(nymph[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.n<-tapply(nymph[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.n<-sd.enhanced.n/sqrt(n.enhanced.n)

#######################################################################################################
#When updating, add an additional label for the new time
#Update the name of the jpeg to the appropriate date


#Line graphs for ACP nymphs
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_nymphs_grove_aug17.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,50), type="n", xlab="Sampling time", ylab="Average ACP nymph count per flush", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July", "Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.n-se.control.n[1:length(unique(time))], 1:length(unique(time)), y.control.n+se.control.n[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.n-se.enhanced.n[1:length(unique(time))], 1:length(unique(time)), y.enhanced.n+se.enhanced.n[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.n, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.n+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1, y=50, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

########################

jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_nymphs_grove_boxplot_aug17.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(grove, nymph, ylim=c(0,50), las=1,xlab="", ylab="", xaxt="n", outline=FALSE, 
     col=c("orange","blue"))
mtext(side=1,line=3, text="Grove", cex=1.5)
mtext(side=2, line=3, text="Avg. ACP nymphs/flush", cex=1.5)
axis(side=1, at=c(1:2),labels=c("Control","Enhanced"))

dev.off()


#######################################################################################################

#Comparing ACP eggs in the enhanced vs. control grove 

y.control.e<-tapply(eggs[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.e<- y.cont[ordered.time]
sd.control.e<-tapply(eggs[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.e<-sd.cont[ordered.time]
n.control.e<-tapply(eggs[which(grove== "control")], time[which(grove== "control")], length)
#n.control.e<-n.cont[ordered.time]
se.control.e<-sd.control.e/sqrt(n.control.e)

y.enhanced.e<-tapply(eggs[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.e<-tapply(eggs[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.e<-tapply(eggs[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.e<-sd.enhanced.e/sqrt(n.enhanced.e)

#######################################################################################################

#Line graphs for ACP eggs
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_eggs_grove_aug17.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,60), type="n", xlab="Sampling time", ylab="Average ACP egg count per flush", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July", "Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.e-se.control.e[1:length(unique(time))], 1:length(unique(time)), y.control.e+se.control.e[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.e-se.enhanced.e[1:length(unique(time))], 1:length(unique(time)), y.enhanced.e+se.enhanced.e[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.e, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.e+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1, y=60, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

#########################

jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_eggs_grove_boxplot_aug17.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(grove, eggs, ylim=c(0,20), las=1,xlab="", ylab="", xaxt="n", outline=FALSE, 
     col=c("orange","blue"))
mtext(side=1,line=3, text="Grove", cex=1.5)
mtext(side=2, line=3, text="Avg. ACP eggs/flush", cex=1.5)
axis(side=1, at=c(1:2),labels=c("Control","Enhanced"))

dev.off()

#######################################################################################################


#Comparing CLM in the enhanced vs. control grove 

y.cont.clm<-tapply(clm[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
y.control.clm<- y.cont.clm[ordered.time]
sd.cont.clm<-tapply(clm[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
sd.control.clm<-sd.cont.clm[ordered.time]
n.cont.clm<-tapply(clm[which(grove== "control")], time[which(grove== "control")], length)
n.control.clm<-n.cont.clm[ordered.time]
se.control.clm<-sd.control.clm/sqrt(n.control.clm)

y.enh.clm<-tapply(clm[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
y.enhanced.clm<- y.enh.clm[ordered.time]
sd.enh.clm<-tapply(clm[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
sd.enhanced.clm<- sd.enh.clm[ordered.time]
n.enh.clm<-tapply(clm[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
n.enhanced.clm<-n.enh.clm[ordered.time]
se.enhanced.clm<-sd.enhanced.clm/sqrt(n.enhanced.clm)

#######################################################################################################

#Line graphs for ACP clm
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_clm_grove_aug17.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,1), type="n", xlab="Sampling time", ylab="Avg. proportion of leaves w/CLM per flush", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July", "Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.clm-se.control.clm[1:length(unique(time))], 1:length(unique(time)), y.control.clm+se.control.clm[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.clm-se.enhanced.clm[1:length(unique(time))], 1:length(unique(time)), y.enhanced.clm+se.enhanced.clm[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.clm, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.clm+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1, y=1, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

#######################################################################################################

jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_clm_grove_boxplot_aug17.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(grove,clm, ylim=c(0,1), las=1,xlab="", ylab="", xaxt="n", outline=FALSE, 
     col=c("orange","blue"))
mtext(side=1,line=3, text="Grove", cex=1.5)
mtext(side=2, line=3, text="Avg. proportion CLM/flush", cex=1.5)
axis(side=1, at=c(1:2),labels=c("Control","Enhanced"))

dev.off()