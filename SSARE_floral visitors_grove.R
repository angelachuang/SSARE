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
vis<-read.csv("C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\SSARE_Floral visitor_Plot level_DATA.csv")
date<-vis$Date
time<-factor(vis$Timepoint, levels = c("April","May","June","July","August"))
grove<-vis$?..Grove
rep<-as.factor(vis$Rep)
plot<-vis$Plot
tx<-vis$Treatment
abun<-vis$Overall_abundance
bee<-vis$All_bees
poll<-vis$Pollinator_abundance
fly<-vis$All_flies

#ordered.time<- factor(c("1","2","3","4", "5")) 

#Comparing average floral visitor abundance per plot in each grove

y.control.abun<-tapply(abun[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.abun<- y.cont[ordered.time]
sd.control.abun<-tapply(abun[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.abun<-sd.cont[ordered.time]
n.control.abun<-tapply(abun[which(grove== "control")], time[which(grove== "control")], length)
#n.control.abun<-n.cont[ordered.time]
se.control.abun<-sd.control.abun/sqrt(n.control.abun)

y.enhanced.abun<-tapply(abun[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.abun<-tapply(abun[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.abun<-tapply(abun[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.abun<-sd.enhanced.abun/sqrt(n.enhanced.abun)

#######################################################################################################

#Line graphs for average visitor abundance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_abundance_grove_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,6), type="n", xlab="Sampling time", ylab="Avg. floral visitor abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July", "August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.abun-se.control.abun[1:length(unique(time))], 1:length(unique(time)), y.control.abun+se.control.abun[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.abun-se.enhanced.abun[1:length(unique(time))], 1:length(unique(time)), y.enhanced.abun+se.enhanced.abun[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.abun, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.abun+0.001, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1.5, y=6, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

#######################################################################################################

#Comparing average bee abundance per plot in each grove

y.control.bee<-tapply(bee[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.bee<- y.cont[ordered.time]
sd.control.bee<-tapply(bee[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.bee<-sd.cont[ordered.time]
n.control.bee<-tapply(bee[which(grove== "control")], time[which(grove== "control")], length)
#n.control.bee<-n.cont[ordered.time]
se.control.bee<-sd.control.bee/sqrt(n.control.bee)

y.enhanced.bee<-tapply(bee[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.bee<-tapply(bee[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.bee<-tapply(bee[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.bee<-sd.enhanced.bee/sqrt(n.enhanced.bee)

#######################################################################################################

#Line graphs for average bee abundance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_bees_grove_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,2), type="n", xlab="Sampling time", ylab="Avg. bee abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July","August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.bee-se.control.bee[1:length(unique(time))], 1:length(unique(time)), y.control.bee+se.control.bee[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.bee-se.enhanced.bee[1:length(unique(time))], 1:length(unique(time)), y.enhanced.bee+se.enhanced.bee[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.bee, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.bee+0.03, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1.5, y=2, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()


#######################################################################################################

#Comparing average fly abundance per plot in each grove

y.control.fly<-tapply(fly[which(grove== "control")], time[which(grove== "control")], mean, na.rm=TRUE)
#y.control.fly<- y.cont[ordered.time]
sd.control.fly<-tapply(fly[which(grove== "control")], time[which(grove== "control")], sd, na.rm=TRUE)
#sd.control.fly<-sd.cont[ordered.time]
n.control.fly<-tapply(fly[which(grove== "control")], time[which(grove== "control")], length)
#n.control.fly<-n.cont[ordered.time]
se.control.fly<-sd.control.fly/sqrt(n.control.fly)

y.enhanced.fly<-tapply(fly[which(grove== "enhanced")], time[which(grove== "enhanced")], mean, na.rm=TRUE)
sd.enhanced.fly<-tapply(fly[which(grove== "enhanced")], time[which(grove== "enhanced")], sd, na.rm=TRUE)
n.enhanced.fly<-tapply(fly[which(grove== "enhanced")], time[which(grove== "enhanced")], length)
se.enhanced.fly<-sd.enhanced.fly/sqrt(n.enhanced.fly)

#######################################################################################################

#Line graphs for average flies at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_flies_grove_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,1.5), type="n", xlab="Sampling time", ylab="Avg. fly abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July","August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control.fly-se.control.fly[1:length(unique(time))], 1:length(unique(time)), y.control.fly+se.control.fly[1:length(unique(time))], col = alpha("orange",0.4), lwd=1)
segments(1:length(unique(time)), y.enhanced.fly-se.enhanced.fly[1:length(unique(time))], 1:length(unique(time)), y.enhanced.fly+se.enhanced.fly[1:length(unique(time))], col = alpha("blue",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.control.fly, type="b", col="orange", bg="orange",lwd=1.5, lty=1, pch=21, cex=2)
lines(c(1:length(unique(time))), y.enhanced.fly+0.03, type="b", col="blue", lwd=1.5, lty=2, pch=18, cex=2)

legend(x=1.5, y=1.5, title="Grove", bty="n", 
       legend =c("Enhanced", "Control"), pt.cex=1.5,
       pch=c(18,21), pt.bg=c("blue","orange"), col=c("blue","orange"))

dev.off()

