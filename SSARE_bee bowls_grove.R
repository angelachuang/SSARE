################################################################################################################
### Project: SSARE Windbreaks ##################################################################################
## Date: July 2021  ## Authors: Chuang, ... ## Year: 2021 ## Version: Meeting with Xavier  #####################
################################################################################################################
### Purpose: Analyze Floral Visitor Data by Grove by Timepoint ##########################################################
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

if (!require(ggpubr)) {  
  install.packages("ggpubr", repos = "http://cran.us.r-project.org")
  require(ggpubr)
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
library(ggpubr) #to make a nice box and whisker plot

###Set up my working directory so the csv ends up in the right folder
setwd("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\")

######################
### Load dataset
beeb<-read.csv("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\SSARE_Bee Bowl_Floral Tx_DATA.csv")
date<-beeb$Date
time<-factor(beeb$Timepoint, levels = c("April","May","June","July", "August"))
grove<-as.factor(beeb$?..Grove)
bees<-as.numeric(beeb$All_bees)
pred<-as.numeric(beeb$All_predators)

#ordered.time<- factor(c("1","2","3","4", "5","6","7","8","9")) 

#Comparing average bee abundance per plot in each grove

y.control.bees<-tapply(bees[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.bees<- y.cont[ordered.time]
sd.control.bees<-tapply(bees[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.bees<-sd.cont[ordered.time]
n.control.bees<-tapply(bees[which(grove== "control")], time[which(grove== "control")], length)
#n.control.bees<-n.cont[ordered.time]
se.control.bees<-sd.control.bees/sqrt(n.control.bees)

y.enhanced.bees<-tapply(bees[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.bees<-tapply(bees[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.bees<-tapply(bees[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.bees<-sd.enhanced.bees/sqrt(n.enhanced.bees)

#######################################################################################################

#Line graphs for bees
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\bee bowls_bees_grove_aug19.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,6), type="n", xlab="Sampling time", ylab="Average bee abundance/sample", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July","August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.bees-se.control.bees[1:length(unique(time))], 1:length(unique(time)), y.control.bees+se.control.bees[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.bees-se.enhanced.bees[1:length(unique(time))], 1:length(unique(time)), y.enhanced.bees+se.enhanced.bees[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.bees, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.bees+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=3, y=6, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

#######################################################################################################

#Bee abundance box plot by grove 

#ggpar(ylim=c(0,5))
#beeplot<-ggboxplot(data=beeb,x="?..Grove", y="All_bees", add="jitter", 
          xlab="Grove", ylab="Bee abundance/plot", binwidth=20)
#ggpar(p=beeplot, ylim=c(0,5))

boxplot(bees~grove, ylim=c(0,5), #notch=TRUE,
        col=(c("blue","orange")), las=1, 
        xlab="Grove Treatment", ylab="Bee abundance/plot", xaxn=c("Control","Enhanced"))
points(x=(rep(1,length(bees[which(grove=="control")]))), y= jitter(bees[which(grove=="control")], factor=1, amount=.05), pch=21)
points(x=(rep(2,length(bees[which(grove=="enhanced")]))), y= jitter(bees[which(grove=="enhanced")], factor=1, amount=.05), pch=21)

plot(grove,bees, ylim=c(0,2), las=1,xlab="Grove treatment", ylab="Bee abundance/plot")

########






y.control.pred<-tapply(pred[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.pred<- y.cont[ordered.time]
sd.control.pred<-tapply(pred[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.pred<-sd.cont[ordered.time]
n.control.pred<-tapply(pred[which(grove== "control")], time[which(grove== "control")], length)
#n.control.pred<-n.cont[ordered.time]
se.control.pred<-sd.control.pred/sqrt(n.control.pred)

y.enhanced.pred<-tapply(pred[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.pred<-tapply(pred[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.pred<-tapply(pred[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.pred<-sd.enhanced.pred/sqrt(n.enhanced.pred)

#######################################################################################################

#Line graphs for predators
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\bee bowls_predators_grove_aug19.jpg", 
     width = 6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,1), type="n", xlab="Sampling time", ylab="Average predator abundance/sample", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July", "August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.pred-se.control.pred[1:length(unique(time))], 1:length(unique(time)), y.control.pred+se.control.pred[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.pred-se.enhanced.pred[1:length(unique(time))], 1:length(unique(time)), y.enhanced.pred+se.enhanced.pred[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.pred, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.pred+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=3, y=1, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()


#####

boxplot(pred~grove, ylim=c(0,5), #notch=TRUE,
        col=(c("blue","orange")), las=1, 
        xlab="Grove Treatment", ylab="Predator abundance/plot", xaxn=c("Control","Enhanced"))
points(x=(rep(1,length(pred[which(grove=="control")]))), y= jitter(pred[which(grove=="control")], factor=1, amount=.05), pch=21)
points(x=(rep(2,length(pred[which(grove=="enhanced")]))), y= jitter(pred[which(grove=="enhanced")], factor=1, amount=.05), pch=21)


plot(grove,pred, ylim=c(0,0.25), las=1,xlab="Grove treatment", ylab="Predator abundance/plot")

