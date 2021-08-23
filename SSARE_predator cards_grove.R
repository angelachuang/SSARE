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

#ordered.time<- factor(c("1","2","3","4", "5","6","7")) 

#Comparing number of intact eggs in the enhanced vs. control grove 

y.control.intact<-tapply(intact[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.intact<- y.cont[ordered.time]
sd.control.intact<-tapply(intact[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.intact<-sd.cont[ordered.time]
n.control.intact<-tapply(intact[which(grove== "control")], time[which(grove== "control")], length)
#n.control.intact<-n.cont[ordered.time]
se.control.intact<-sd.control.intact/sqrt(n.control.intact)

y.enhanced.intact<-tapply(intact[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.intact<-tapply(intact[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.intact<-tapply(intact[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.intact<-sd.enhanced.intact/sqrt(n.enhanced.intact)

#######################################################################################################

#Line graphs for intact eggs from predator cards
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\predator cards_intact_grove_jul21.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,1), type="n", xlab="Sampling time", ylab="Avg proportion of intact eggs", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.intact-se.control.intact[1:length(unique(time))], 1:length(unique(time)), y.control.intact+se.control.intact[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.intact-se.enhanced.intact[1:length(unique(time))], 1:length(unique(time)), y.enhanced.intact+se.enhanced.intact[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.intact, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.intact+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=3, y=1, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()