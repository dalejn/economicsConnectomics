rm(list=ls())

################################
# Load workspace and libraries #
################################

load('0_finalData.RData')

library(ggplot2)
library(ppcor)
library(mgcv)
library(visreg)
library(scatterplot3d)

#####################################
# Make example rate-distortion plot #
#####################################

par(xaxs="i", yaxs="i") 
d <- seq(from=1e-3,to=1,length.out=100)
r <- (1/2) * log2(1/d)
plot(c(0,1),c(0,4),type="n",
     xlab="Distortion",
     ylab="Rate (bits)",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(d,r,col="red",lwd=3)


setEPS()
postscript("figures/fig4/toyRD.eps")
par(xaxs="i", yaxs="i") 
plot(c(0,1),c(0,4),type="n",
     xlab="Distortion",
     ylab="Rate (bits)",
     cex.lab=1.5, cex.axis=1.5, 
     cex.main=1.5, cex.sub=1.5)
lines(d,r,col="red",lwd=3)

dev.off()

##############################################################
# Plot brain network resource efficiency with random network #
##############################################################

# see MATLAB code fig4_randomVsBrainRD.m

################################################
# Plot compression efficiency with development #
################################################


unbiasedSlopes <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/unbiasedSlopes.txt', sep=' ',header=F))
colnames(unbiasedSlopes) <- c('scanid','unbiasedSlopes')
QA_df <- merge(QA_df, unbiasedSlopes, by=c('scanid'))


QA_df$ageInYears <- QA_df$ageAtScan1/12
globEffCbf<- gam(unbiasedSlopes ~ s(ageInYears,k=4) + sex+s(ageInYears, by=sex, k=4)+degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df)
summary(globEffCbf)
levels(QA_df$sex)<- c("2","1")

setEPS()
postscript("figures/fig4/channelCapacityByAge.eps")
visreg(globEffCbf,"ageInYears","sex",overlay=T,xlim=c(5,25),gg=T)+theme_classic(base_size=20)
dev.off()