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

###############################################################
# Individual fits for RD of unbiased random walk & chemotaxis #
###############################################################

# See MATLAB code fig5_fitRD.m

#############################################
# Compare unbiased random walk & chemotaxis #
#############################################

for (i in c(1,2,3,4,5,6,7,8,9,92,94,96,98,999))
{
  filepath <- paste0('data/glasser__resource_efficiency0',i,'.txt')
  resEff <- read.csv(filepath,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
  is.na(resEff) <- do.call(cbind,lapply(resEff, is.infinite))
  resEff <- cbind(resEff[1],rowMeans(resEff[2:ncol(resEff)],na.rm=TRUE))
  resEff <- resEff[which(resEff[,2]>0),]
  if (i==999)
  {
    resEff$probability <- rep(i/1000,length(resEff[,1]))
  } else if (i>90 & i <100) {
    resEff$probability <- rep(i/100,length(resEff[,1]))
  } else {
    resEff$probability <- rep(i/10,length(resEff[,1]))
  }
  resEff$cbfBiased <- rep(0,length(resEff[,1]))
  colnames(resEff) <- c('scanid',"resources","probability","cbfBiased")
  assign(paste0("df",i),as.data.frame(resEff))
  
  filepath2 <- paste0('/data/glasser_cbfBiasedAttractive_resource_efficiency0',i,'.txt')
  resEff <- read.csv(filepath2,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
  is.na(resEff) <- do.call(cbind,lapply(resEff, is.infinite))
  resEff <- cbind(resEff[1],rowMeans(resEff[2:ncol(resEff)],na.rm=TRUE))
  resEff <- resEff[which(resEff[,2]>0),]
  if (i==999)
  {
    resEff$probability <- rep(i/1000,length(resEff[,1]))
  } else if (i>90 & i <100) {
    resEff$probability <- rep(i/100,length(resEff[,1]))
  } else {
    resEff$probability <- rep(i/10,length(resEff[,1]))
  }
  resEff$cbfBiased <- rep(1,length(resEff[,1]))
  colnames(resEff) <- c('scanid',"resources","probability","cbfBiased")
  assign(paste0("df",i,"biased"),as.data.frame(resEff))
  
  filepath3 <- paste0('/data/glasser_cbfBiasedRepulsive_resource_efficiency0',i,'.txt')
  resEff <- read.csv(filepath3,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
  is.na(resEff) <- do.call(cbind,lapply(resEff, is.infinite))
  resEff <- cbind(resEff[1],rowMeans(resEff[2:ncol(resEff)],na.rm=TRUE))
  resEff <- resEff[which(resEff[,2]>0),]
  if (i==999)
  {
    resEff$probability <- rep(i/1000,length(resEff[,1]))
  } else if (i>90 & i <100) {
    resEff$probability <- rep(i/100,length(resEff[,1]))
  } else {
    resEff$probability <- rep(i/10,length(resEff[,1]))
  }
  resEff$cbfBiased <- rep(2,length(resEff[,1]))
  colnames(resEff) <- c('scanid',"resources","probability","cbfBiased")
  assign(paste0("df",i,"biasedRepulse"),as.data.frame(resEff))
}

df1 <- rbind(df1,df1biased,df1biasedRepulse,df2,df2biased,df2biasedRepulse,df3,df3biased,df3biasedRepulse
             ,df4,df4biased,df4biasedRepulse,df5,df5biased,df5biasedRepulse,df6,df6biased,df6biasedRepulse,df7,df7biased,df7biasedRepulse,
             df8,df8biased,df8biasedRepulse,df9,df9biased,df9biasedRepulse,df92,df92biased,df92biasedRepulse,df94,df94biased,df94biasedRepulse,
             df96,df96biased,df96biasedRepulse,df98,df98biased,df98biasedRepulse,df999,df999biased,df999biasedRepulse)


df1$cbfBiased <- as.factor(df1$cbfBiased+1)
df1$probability <- as.factor(df1$probability)

p <- ggplot(df1,aes(probability,resources,fill=cbfBiased))+geom_boxplot(outlier.size=0.05)+
  facet_wrap(~probability, scale="free",ncol=3)+
  theme_classic(base_size=20)
p <- p+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
p

ggsave('figures/fig5/boxWhiskerPlots.eps',device='eps',width=7.18,height=6.31)

# Statistical testing

summary(lm(resources~ probability+cbfBiased+ probability*cbfBiased, data=df1))

df2 <- df1[which(df1$probability=="0.999"),]
t.test(df2[which(df2$cbfBiased==1),][,2], df2[which(df2$cbfBiased==2),][,2])
t.test(df2[which(df2$cbfBiased==1),][,2], df2[which(df2$cbfBiased==3),][,2])
t.test(df2[which(df2$cbfBiased==2),][,2], df2[which(df2$cbfBiased==3),][,2])