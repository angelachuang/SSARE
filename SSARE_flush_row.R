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
row<-as.factor(flush$Row)
tree<-flush$Tree
tree2<-flush$T
set<-flush$Set
tx<-flush$Treatment
nymph<-flush$ACP_nymph
eggs<-flush$ACP_egg
clm<-as.numeric(flush$prop_CLM)

ordered.time<- factor(c("1","2","3","4", "5","6","7","8","9")) 

#######################################################################################################

#Comparing eggs per row

y.1.eggs<-tapply(eggs[which(row=="1")], time[which(row=="1")], mean, na.rm=TRUE)
sd.1.eggs<-tapply(eggs[which(row=="1")], time[which(row=="1")], sd, na.rm=TRUE)
n.1.eggs<-tapply(eggs[which(row=="1")], time[which(row=="1")], length)
se.1.eggs<-sd.1.eggs/sqrt(n.1.eggs)

y.3.eggs<-tapply(eggs[which(row=="3")], time[which(row=="3")], mean, na.rm=TRUE)
sd.3.eggs<-tapply(eggs[which(row=="3")], time[which(row=="3")], sd, na.rm=TRUE)
n.3.eggs<-tapply(eggs[which(row=="3")], time[which(row=="3")], length)
se.3.eggs<-sd.3.eggs/sqrt(n.3.eggs)

y.6.eggs<-tapply(eggs[which(row=="6")], time[which(row=="6")], mean, na.rm=TRUE)
sd.6.eggs<-tapply(eggs[which(row=="6")], time[which(row=="6")], sd, na.rm=TRUE)
n.6.eggs<-tapply(eggs[which(row=="6")], time[which(row=="6")], length)
se.6.eggs<-sd.6.eggs/sqrt(n.6.eggs)


#######################################################################################################

#Line graphs for average visitor eggsdance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_eggs_row_aug17.jpg", 
     width =6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,80), type="n", xlab="Sampling time", ylab="Avg. eggs abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.1.eggs-se.1.eggs[1:length(unique(time))], 1:length(unique(time)), y.1.eggs+se.1.eggs[1:length(unique(time))], col = alpha("firebrick2",0.4) , lwd=1)
segments(1:length(unique(time)), y.3.eggs-se.3.eggs[1:length(unique(time))], 1:length(unique(time)), y.3.eggs+se.3.eggs[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.6.eggs-se.6.eggs[1:length(unique(time))], 1:length(unique(time)), y.6.eggs+se.6.eggs[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.1.eggs, type="b", col="firebrick2", bg="firebrick2", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.3.eggs, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.6.eggs, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=80, title="Distance from windbreak", bty="n", 
       legend =c("Row 1","Row 3","Row 6"), pt.cex=1.5,
       pch=c(21), pt.bg=c("firebrick2","darkgoldenrod2","aquamarine4"), 
       col=c("firebrick2","darkgoldenrod2","aquamarine4"))

dev.off()


#######################################################################################################

#Comparing nymph per row

y.1.nymph<-tapply(nymph[which(row=="1")], time[which(row=="1")], mean, na.rm=TRUE)
sd.1.nymph<-tapply(nymph[which(row=="1")], time[which(row=="1")], sd, na.rm=TRUE)
n.1.nymph<-tapply(nymph[which(row=="1")], time[which(row=="1")], length)
se.1.nymph<-sd.1.nymph/sqrt(n.1.nymph)

y.3.nymph<-tapply(nymph[which(row=="3")], time[which(row=="3")], mean, na.rm=TRUE)
sd.3.nymph<-tapply(nymph[which(row=="3")], time[which(row=="3")], sd, na.rm=TRUE)
n.3.nymph<-tapply(nymph[which(row=="3")], time[which(row=="3")], length)
se.3.nymph<-sd.3.nymph/sqrt(n.3.nymph)

y.6.nymph<-tapply(nymph[which(row=="6")], time[which(row=="6")], mean, na.rm=TRUE)
sd.6.nymph<-tapply(nymph[which(row=="6")], time[which(row=="6")], sd, na.rm=TRUE)
n.6.nymph<-tapply(nymph[which(row=="6")], time[which(row=="6")], length)
se.6.nymph<-sd.6.nymph/sqrt(n.6.nymph)

#######################################################################################################

#Line graphs for nymphs at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_nymph_row_aug17.jpg", 
     width =6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,50), type="n", xlab="Sampling time", ylab="Avg. nymph abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.1.nymph-se.1.nymph[1:length(unique(time))], 1:length(unique(time)), y.1.nymph+se.1.nymph[1:length(unique(time))], col = alpha("firebrick2",0.4) , lwd=1)
segments(1:length(unique(time)), y.3.nymph-se.3.nymph[1:length(unique(time))], 1:length(unique(time)), y.3.nymph+se.3.nymph[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.6.nymph-se.6.nymph[1:length(unique(time))], 1:length(unique(time)), y.6.nymph+se.6.nymph[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.1.nymph, type="b", col="firebrick2", bg="firebrick2", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.3.nymph, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.6.nymph, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=50, title="Distance from windbreak", bty="n", 
       legend =c("Row 1","Row 3","Row 6"), pt.cex=1.5,
       pch=c(21), pt.bg=c("firebrick2","darkgoldenrod2","aquamarine4"), 
       col=c("firebrick2","darkgoldenrod2","aquamarine4"))

dev.off()



#######################################################################################################

#Line graphs for average visitor eggsdance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_eggs_row_aug17.jpg", 
     width =6, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,80), type="n", xlab="Sampling time", ylab="Avg. eggs abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July","Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.1.eggs-se.1.eggs[1:length(unique(time))], 1:length(unique(time)), y.1.eggs+se.1.eggs[1:length(unique(time))], col = alpha("firebrick2",0.4) , lwd=1)
segments(1:length(unique(time)), y.3.eggs-se.3.eggs[1:length(unique(time))], 1:length(unique(time)), y.3.eggs+se.3.eggs[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.6.eggs-se.6.eggs[1:length(unique(time))], 1:length(unique(time)), y.6.eggs+se.6.eggs[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.1.eggs, type="b", col="firebrick2", bg="firebrick2", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.3.eggs, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.6.eggs, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=80, title="Distance from windbreak", bty="n", 
       legend =c("Row 1","Row 3","Row 6"), pt.cex=1.5,
       pch=c(21), pt.bg=c("firebrick2","darkgoldenrod2","aquamarine4"), 
       col=c("firebrick2","darkgoldenrod2","aquamarine4"))

dev.off()


#######################################################################################################

#Comparing clm per row

y.1.clm<-tapply(clm[which(row=="1")], time[which(row=="1")], mean, na.rm=TRUE)
sd.1.clm<-tapply(clm[which(row=="1")], time[which(row=="1")], sd, na.rm=TRUE)
n.1.clm<-tapply(clm[which(row=="1")], time[which(row=="1")], length)
se.1.clm<-sd.1.clm/sqrt(n.1.clm)

y.3.clm<-tapply(clm[which(row=="3")], time[which(row=="3")], mean, na.rm=TRUE)
sd.3.clm<-tapply(clm[which(row=="3")], time[which(row=="3")], sd, na.rm=TRUE)
n.3.clm<-tapply(clm[which(row=="3")], time[which(row=="3")], length)
se.3.clm<-sd.3.clm/sqrt(n.3.clm)

y.6.clm<-tapply(clm[which(row=="6")], time[which(row=="6")], mean, na.rm=TRUE)
sd.6.clm<-tapply(clm[which(row=="6")], time[which(row=="6")], sd, na.rm=TRUE)
n.6.clm<-tapply(clm[which(row=="6")], time[which(row=="6")], length)
se.6.clm<-sd.6.clm/sqrt(n.6.clm)

#######################################################################################################

#Line graphs for average visitor clmdance at the plot level for each grove
jpeg(filename="C:\\Users\\chuan\\Dropbox\\Ongoing research\\Diepenbrock Lab\\SSARE Windbreaks\\Data\\Data visualizations\\flush_clm_row_aug17.jpg", 
     width =12, height = 5.5, units= "in" ,res=300)

par(mar=c(5,5,1,2))
plot(1, xlim= c(1,length(unique(time))), ylim = c(0,.8), type="n", xlab="Sampling time", ylab="Avg. CLM abundance per plot", xaxt="n", las=1, cex.lab=1.5)
axis(side=1, at= 1:length(unique(time)), labels=c("Mid-April" , "Late April" , "Mid-May" , "Late May" , "Mid-June", "Late June","Mid-July", "Late July","Mid-August"), cex.lab=1.5)

segments(1:length(unique(time)), y.1.clm-se.1.clm[1:length(unique(time))], 1:length(unique(time)), y.1.clm+se.1.clm[1:length(unique(time))], col = alpha("firebrick2",0.4) , lwd=1)
segments(1:length(unique(time)), y.3.clm-se.3.clm[1:length(unique(time))], 1:length(unique(time)), y.3.clm+se.3.clm[1:length(unique(time))], col = alpha("darkgoldenrod2",0.4) , lwd=1)
segments(1:length(unique(time)), y.6.clm-se.6.clm[1:length(unique(time))], 1:length(unique(time)), y.6.clm+se.6.clm[1:length(unique(time))], col = alpha("aquamarine4",0.4) , lwd=1)

lines(c(1:length(unique(time))), y.1.clm, type="b", col="firebrick2", bg="firebrick2", lwd=1.5, lty=2, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.3.clm, type="b", col="darkgoldenrod2", bg="darkgoldenrod2",lwd=1.5, lty=1, pch=21, cex=1.5)
lines(c(1:length(unique(time))), y.6.clm, type="b", col="aquamarine4", bg="aquamarine4", lwd=1.5, lty=2, pch=21, cex=1.5)


legend(x=1, y=.8, title="Distance from windbreak", bty="n", 
       legend =c("Row 1","Row 3","Row 6"), pt.cex=1.5,
       pch=c(21), pt.bg=c("firebrick2","darkgoldenrod2","aquamarine4"), 
       col=c("firebrick2","darkgoldenrod2","aquamarine4"))

dev.off()