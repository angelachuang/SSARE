################################################################################################################
### Project: SSARE Windbreaks ##################################################################################
## Date: July 2021  ## Authors: Chuang, ... ## Year: 2021 ## Version: Meeting with Xavier  #####################
################################################################################################################
### Purpose: Analyze Floral Visitor Data by Treatment by Timepoint ##########################################################
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
tx<-factor(vis$Treatment, levels=c("control_grove","control","bush","vine","flower","all"))
abun<-as.numeric(vis$Overall_abundance)
bee<-as.numeric(vis$All_bees)
poll<-as.numeric(vis$Pollinator_abundance)
fly<-as.numeric(vis$All_flies)

#ordered.time<- factor(c("1","2","3","4", "5","6","7")) 

#######################################################################################################


#Comparing average floral visitor abundance per treatment

y.control_grove.abun<-tapply(abun[which(tx=="control_grove")], time[which(tx=="control_grove")], mean, na.rm=TRUE)
#y.control_grove.abun<- y.cont[ordered.time]
sd.control_grove.abun<-tapply(abun[which(tx=="control_grove")], time[which(tx=="control_grove")], sd, na.rm=TRUE)
#sd.control_grove.abun<-sd.cont[ordered.time]
n.control_grove.abun<-tapply(abun[which(tx=="control_grove")], time[which(tx=="control_grove")], length)
#n.control_grove.abun<-n.cont[ordered.time]
se.control_grove.abun<-sd.control_grove.abun/sqrt(n.control_grove.abun)

y.control.abun<-tapply(abun[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.abun<-tapply(abun[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.abun<-tapply(abun[which(tx=="control")], time[which(tx=="control")], length)
se.control.abun<-sd.control.abun/sqrt(n.control.abun)

y.flower.abun<-tapply(abun[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.abun<-tapply(abun[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.abun<-tapply(abun[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.abun<-sd.flower.abun/sqrt(n.flower.abun)

y.vine.abun<-tapply(abun[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.abun<-tapply(abun[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.abun<-tapply(abun[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.abun<-sd.vine.abun/sqrt(n.vine.abun)

y.bush.abun<-tapply(abun[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.abun<-tapply(abun[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.abun<-tapply(abun[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.abun<-sd.bush.abun/sqrt(n.bush.abun)

y.all.abun<-tapply(abun[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.abun<-tapply(abun[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.abun<-tapply(abun[which(tx=="all")], time[which(tx=="all")], length)
se.all.abun<-sd.all.abun/sqrt(n.all.abun)

#######################################################################################################

#Line graphs for average visitor abundance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_abundance_tx_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,15), type="n", xlab="Sampling time", ylab="Avg. floral visitor abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July", "August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control_grove.abun-se.control_grove.abun[1:length(unique(time))], 1:length(unique(time)), y.control_grove.abun+se.control_grove.abun[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.abun-se.control.abun[1:length(unique(time))], 1:length(unique(time)), y.control.abun+se.control.abun[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.abun-se.flower.abun[1:length(unique(time))], 1:length(unique(time)), y.flower.abun+se.flower.abun[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.abun-se.vine.abun[1:length(unique(time))], 1:length(unique(time)), y.vine.abun+se.vine.abun[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.abun-se.bush.abun[1:length(unique(time))], 1:length(unique(time)), y.bush.abun+se.bush.abun[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.abun-se.all.abun[1:length(unique(time))], 1:length(unique(time)), y.all.abun+se.all.abun[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control_grove.abun, type="b", col="indianred4", bg="indianred4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.abun, type="b", col="pink", bg="pink", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.abun, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.abun, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.abun, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.abun, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=2, y=15, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"), 
       col=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"))

dev.off()

####

plot(tx,abun, ylim=c(0,20), las=1,xlab="Grove treatment", ylab="Floral visitor abundance/plot")


#######################################################################################################

#Comparing bees per treatment

y.control_grove.bee<-tapply(bee[which(tx=="control_grove")], time[which(tx=="control_grove")], mean, na.rm=TRUE)
#y.control_grove.bee<- y.cont[ordered.time]
sd.control_grove.bee<-tapply(bee[which(tx=="control_grove")], time[which(tx=="control_grove")], sd, na.rm=TRUE)
#sd.control_grove.bee<-sd.cont[ordered.time]
n.control_grove.bee<-tapply(bee[which(tx=="control_grove")], time[which(tx=="control_grove")], length)
#n.control_grove.bee<-n.cont[ordered.time]
se.control_grove.bee<-sd.control_grove.bee/sqrt(n.control_grove.bee)

y.control.bee<-tapply(bee[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.bee<-tapply(bee[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.bee<-tapply(bee[which(tx=="control")], time[which(tx=="control")], length)
se.control.bee<-sd.control.bee/sqrt(n.control.bee)

y.flower.bee<-tapply(bee[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.bee<-tapply(bee[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.bee<-tapply(bee[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.bee<-sd.flower.bee/sqrt(n.flower.bee)

y.vine.bee<-tapply(bee[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.bee<-tapply(bee[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.bee<-tapply(bee[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.bee<-sd.vine.bee/sqrt(n.vine.bee)

y.bush.bee<-tapply(bee[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.bee<-tapply(bee[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.bee<-tapply(bee[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.bee<-sd.bush.bee/sqrt(n.bush.bee)

y.all.bee<-tapply(bee[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.bee<-tapply(bee[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.bee<-tapply(bee[which(tx=="all")], time[which(tx=="all")], length)
se.all.bee<-sd.all.bee/sqrt(n.all.bee)

#######################################################################################################

#Line graphs for average visitor beedance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_bees_tx_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,4), type="n", xlab="Sampling time", ylab="Avg. bee abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July", "August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control_grove.bee-se.control_grove.bee[1:length(unique(time))], 1:length(unique(time)), y.control_grove.bee+se.control_grove.bee[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.bee-se.control.bee[1:length(unique(time))], 1:length(unique(time)), y.control.bee+se.control.bee[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.bee-se.flower.bee[1:length(unique(time))], 1:length(unique(time)), y.flower.bee+se.flower.bee[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.bee-se.vine.bee[1:length(unique(time))], 1:length(unique(time)), y.vine.bee+se.vine.bee[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.bee-se.bush.bee[1:length(unique(time))], 1:length(unique(time)), y.bush.bee+se.bush.bee[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.bee-se.all.bee[1:length(unique(time))], 1:length(unique(time)), y.all.bee+se.all.bee[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control_grove.bee, type="b", col="black", bg="black",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.bee, type="b", col="grey", bg="grey", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.bee, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.bee, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.bee-0.025, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.bee-0.05, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=2, y=4, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"), 
       col=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"))

dev.off()

#########################
#Boxplot of treatments x bees
#Outline = false takes out the outliers. I've done this to make the plots simpler 

jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors boxplot_bees_tx_aug19.jpg", 
     width = 8, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(tx,bee, ylim=c(0,6), las=1,xlab="", ylab="", xaxt="n", outline=FALSE, 
     col=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"))
mtext(side=1,line=3, text="Grove treatment", cex=1.5)
mtext(side=2, line=3, text="Average bee abundance/plot", cex=1.5)
axis(side=1, at=c(1:6),labels=c("Control Grove","Control","Bush","Vine","Wildflower","All"))
abline(v=1.5, lwd=2, lty=2, col="grey40")

dev.off()

#######################################################################################################

#Comparing flies per treatment

y.control_grove.fly<-tapply(fly[which(tx=="control_grove")], time[which(tx=="control_grove")], mean, na.rm=TRUE)
#y.control_grove.fly<- y.cont[ordered.time]
sd.control_grove.fly<-tapply(fly[which(tx=="control_grove")], time[which(tx=="control_grove")], sd, na.rm=TRUE)
#sd.control_grove.fly<-sd.cont[ordered.time]
n.control_grove.fly<-tapply(fly[which(tx=="control_grove")], time[which(tx=="control_grove")], length)
#n.control_grove.fly<-n.cont[ordered.time]
se.control_grove.fly<-sd.control_grove.fly/sqrt(n.control_grove.fly)

y.control.fly<-tapply(fly[which(tx=="control")], time[which(tx=="control")], mean, na.rm=TRUE)
sd.control.fly<-tapply(fly[which(tx=="control")], time[which(tx=="control")], sd, na.rm=TRUE)
n.control.fly<-tapply(fly[which(tx=="control")], time[which(tx=="control")], length)
se.control.fly<-sd.control.fly/sqrt(n.control.fly)

y.flower.fly<-tapply(fly[which(tx=="flower")], time[which(tx=="flower")], mean, na.rm=TRUE)
sd.flower.fly<-tapply(fly[which(tx=="flower")], time[which(tx=="flower")], sd, na.rm=TRUE)
n.flower.fly<-tapply(fly[which(tx=="flower")], time[which(tx=="flower")], length)
se.flower.fly<-sd.flower.fly/sqrt(n.flower.fly)

y.vine.fly<-tapply(fly[which(tx=="vine")], time[which(tx=="vine")], mean, na.rm=TRUE)
sd.vine.fly<-tapply(fly[which(tx=="vine")], time[which(tx=="vine")], sd, na.rm=TRUE)
n.vine.fly<-tapply(fly[which(tx=="vine")], time[which(tx=="vine")], length)
se.vine.fly<-sd.vine.fly/sqrt(n.vine.fly)

y.bush.fly<-tapply(fly[which(tx=="bush")], time[which(tx=="bush")], mean, na.rm=TRUE)
sd.bush.fly<-tapply(fly[which(tx=="bush")], time[which(tx=="bush")], sd, na.rm=TRUE)
n.bush.fly<-tapply(fly[which(tx=="bush")], time[which(tx=="bush")], length)
se.bush.fly<-sd.bush.fly/sqrt(n.bush.fly)

y.all.fly<-tapply(fly[which(tx=="all")], time[which(tx=="all")], mean, na.rm=TRUE)
sd.all.fly<-tapply(fly[which(tx=="all")], time[which(tx=="all")], sd, na.rm=TRUE)
n.all.fly<-tapply(fly[which(tx=="all")], time[which(tx=="all")], length)
se.all.fly<-sd.all.fly/sqrt(n.all.fly)

#######################################################################################################

#Line graphs for average fly abundance at the plot level 
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors_flies_tx_aug19.jpg", 
     width = 4, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,3), type="n", xlab="Sampling time", ylab="Avg. fly abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("April","May","June","July","August"), cex.lab=1.5)

segments(1:length(unique(time)), y.control_grove.fly-se.control_grove.fly[1:length(unique(time))], 1:length(unique(time)), y.control_grove.fly+se.control_grove.fly[1:length(unique(time))], col = alpha("black",0.4), lwd=1)
segments(1:length(unique(time)), y.control.fly-se.control.fly[1:length(unique(time))], 1:length(unique(time)), y.control.fly+se.control.fly[1:length(unique(time))], col = alpha("grey",0.4) , lwd=1)
segments(1:length(unique(time)), y.flower.fly-se.flower.fly[1:length(unique(time))], 1:length(unique(time)), y.flower.fly+se.flower.fly[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.vine.fly-se.vine.fly[1:length(unique(time))], 1:length(unique(time)), y.vine.fly+se.vine.fly[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)
segments(1:length(unique(time)), y.bush.fly-se.bush.fly[1:length(unique(time))], 1:length(unique(time)), y.bush.fly+se.bush.fly[1:length(unique(time))], col = alpha("lightsalmon4",0.4) , lwd=1)
segments(1:length(unique(time)), y.all.fly-se.all.fly[1:length(unique(time))], 1:length(unique(time)), y.all.fly+se.all.fly[1:length(unique(time))], col = alpha("dodgerblue1",0.4) , lwd=1)


lines(c(1:length(unique(time))), y.control_grove.fly, type="b", col="black", bg="black",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.control.fly, type="b", col="grey", bg="grey", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.flower.fly, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.vine.fly, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.bush.fly-0.025, type="b", col="lightsalmon4", bg="lightsalmon4",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.all.fly-0.05, type="b", col="dodgerblue1", bg="dodgerblue1", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=2, y=3, title="Treatment", bty="n", 
       legend =c("Control grove","Control","Wildflower","Vine","Bush","All"), pt.cex=1.5,
       pch=c(21), pt.bg=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"), 
       col=c("black","grey","darkgoldenrod2","aquamarine4","lightsalmon4","dodgerblue1"))

dev.off()


#########################
#Boxplot of treatments x flies
#Outline = false takes out the outliers. I've done this to make the plots simpler 

jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\floral visitors boxplot_flies_tx_aug19.jpg", 
     width = 8, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(tx,fly, ylim=c(0,4), las=1,xlab="", ylab="", xaxt="n", outline=FALSE, 
     col=c("indianred4", "pink","orange2","goldenrod","aquamarine4","royalblue2"))
mtext(side=1,line=3, text="Grove treatment", cex=1.5)
mtext(side=2, line=3, text="Average fly abundance/plot", cex=1.5)
axis(side=1, at=c(1:6),labels=c("Control Grove","Control","Bush","Vine","Wildflower","All"))
abline(v=1.5, lwd=2, lty=2, col="grey40")

dev.off()
