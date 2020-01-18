rm(list=ls())

################################
# Load workspace and libraries #
################################

load('data/0_finalData.RData')
library(ggplot2)
library(ppcor)
library(mgcv)
library(visreg)
library(scatterplot3d)

pathtrans<- read.csv("/data/pathTransArray.txt",sep=" ",header=F)
pathtrans <- cbind(pathtrans[1], rowMeans(pathtrans[2:361], na.rm=T))
colnames(pathtrans)<- c('scanid','pathtrans')

QA_df <- merge(QA_df,pathtrans,by=c("scanid"))

############################################################
# Trade-offs between modularity	and diffusion architecture #
############################################################

# Path transitivity and CBF
globEffCbf<- gam(globalCBF ~ pathtrans + s(ageAtScan1,k=4)+ degree + density + s(ageAtScan1,by=sex,k=4) + sex + meanMotion,fx=TRUE,method="REML",data=QA_df)
visreg(globEffCbf,"pathtrans", xlab="Path transitivity", ylab="Global CBF", gg=T)+theme_classic(base_size=20)+labs(y="Global CBF (ml/100g/min)", x="Path transitivity")
ggsave('figures/fig3/globalCbfPathTransitivity.eps',device='eps',width=7.18,height=6.31)

# Path transitivity and CBF moderated by age
globEffCbf<- gam(globalCBF ~ s(pathtrans,by=ageAtScan1) + s(ageAtScan1,k=4)+ degree + density + s(ageAtScan1,by=sex,k=4) + sex + meanMotion,fx=TRUE,method="REML",data=QA_df)

# Modularity interaction
globEffCbf<- gam(globalCBF ~ pathtrans + Q + pathtrans*Q + s(ageAtScan1,k=4)+ degree + density + s(ageAtScan1,by=sex,k=4) + sex + meanMotion,fx=TRUE,method="REML",data=QA_df)

# Modularity and path transitivity fitness landscape
setEPS()
postscript("figures/fig3/globalCbfPathTransitivityModularity_surface.eps")
vis.gam(globEffCbf,view=c("Q","pathtrans"),type="response",color="topo",theta=235,phi=5, ticktype="detailed",xlab="Modularity",ylab="Path transitivity",zlab="Global CBF")
dev.off()

# Histograms
ggplot(QA_df,aes(Q)) + geom_histogram(bins=20)+theme_classic(base_size=30)+labs(x="Modularity",y="Frequency")+geom_vline(xintercept=mean(QA_df$Q),size=3,linetype="longdash",color="blue")
ggsave('figures/fig3/modularityHistogram.eps',device='eps',width=7.18,height=6.31)

ggplot(QA_df,aes(pathtrans)) + geom_histogram(bins=20)+theme_classic(base_size=30)+labs(x="Path transitivity",y="Frequency")+geom_vline(xintercept=mean(QA_df$pathtrans),size=3,linetype="longdash",color="blue")
ggsave('figures/fig3/pathtransHistogram.eps',device='eps',width=7.18,height=6.31)