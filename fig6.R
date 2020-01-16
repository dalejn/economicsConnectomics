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
library(scales)

####################################
# Lausanne atlas for HIGH fidelity #
####################################

# all individual data points
for (j in c(83,129,234,463,1015))
{
  for (i in c(999,9,5))
  {
    filepath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/lausanneResourceEfficiency/largestComponent/Lausanne',j,'_resource_efficiency0',i,'.txt')
    
    degreePath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/degreeLausanne',j,'.txt')
    degree <- read.csv(degreePath,header=F,sep=' ')
    degree_name <- paste0('degree',i,'_',j,'nodes')
    colnames(degree) <- c("scanid",degree_name)
    
    densityPath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/densityLausanne',j,'.txt')
    density <- read.csv(degreePath,header=F,sep=' ')
    density_name <- paste0('density',i,'_',j,'nodes')
    colnames(density) <- c("scanid",density_name)
    
    resEff <- read.csv(filepath,header=F,sep=' ')
    resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
    is.na(resEff) <- do.call(cbind,lapply(resEff, is.infinite))
    resEff <- cbind(resEff[1],rowMeans(resEff[2:ncol(resEff)],na.rm=TRUE))
    column_new <- paste0('resEff',i,'_',j,'nodes')
    colnames(resEff) <- c('scanid',column_new)
    resEff <- merge(resEff, degree, by=c("scanid"))
    resEff <- merge(resEff, density,by=c("scanid"))
    
    QA_df<- merge(QA_df,resEff,by=c('scanid'))
    QA_minusOutlier <- QA_df[!(QA_df[[column_new]] < mean(QA_df[[column_new]]) - sd(QA_df[[column_new]])*3 | QA_df[[column_new]] > mean(QA_df[[column_new]])+ sd(QA_df[[column_new]])*3),]
    placeholder <- QA_minusOutlier[[column_new]]
    placeholder2 <- QA_minusOutlier[[degree_name]]
    placeholder3 <- QA_minusOutlier[[density_name]]
    placeholder4 <- QA_minusOutlier$ageAtScan1
    placeholder5 <- QA_minusOutlier$sex
    placeholder6 <- QA_minusOutlier$dti64MeanRelRMS
    if (i==999)
    {
      df1 <- as.data.frame(cbind(rep(j,length(placeholder)),rep(0.1,length(placeholder)),placeholder,placeholder2,placeholder3, placeholder4,placeholder5,placeholder6))
    } else{
      df1 <- as.data.frame(cbind(rep(j,length(placeholder)),rep(100-(as.numeric(i)*10),length(placeholder)),placeholder,placeholder2,placeholder3, placeholder4,placeholder5,placeholder6))
    }
    
    if (j==83 & i==999)
    {
      df2 <- df1
    } else{
      df2 <- rbind(df2,df1)
    }
  }
}

df2 <- as.data.frame(df2)
colnames(df2) <- c('nodes','distortion','resources','degree','density','age','sex','motion')
df2$distortion <- as.factor(df2$distortion)
df2$sex <- as.factor(df2$sex)

plotDf <- df2
lm1 <- gam(resources~nodes+degree+density+age+sex+age*sex+motion,method="REML",data=df2)

require(scales)
ggplot(df2,aes(nodes,resources,col=distortion))+geom_point(shape=21,size=.1)+geom_jitter(width=.5)+geom_smooth(method='lm')+theme_classic(base_size=20)
ggsave('figures/fig6/highFidelityRegime.eps',device='eps',width=7.18,height=6.31)

##################################
# Glasser atlas for LOW fidelity #
##################################

# load path complexity 

df1 <- read.csv('/data/jux/BBL/projects/ASLnetwork/results/pathSize.txt',sep=' ',header=F)
colnames(df1)<- c('scanid','pathComplexity')
QA_df <- merge(QA_df, df1,by='scanid')

# Res Eff Glasser

for (i in c(1,2,3,4,5,6,7,8,9,92,94,96,98,999))
{
  filepath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/glasser_resource_efficiency0',i,'.txt')
  resEff <- read.csv(filepath,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
  is.na(resEff) <- do.call(cbind,lapply(resEff, is.infinite))
  resEff <- cbind(resEff[1],rowMeans(resEff[2:ncol(resEff)],na.rm=TRUE))
  colnames(resEff) <- c('scanid',"resources")
  resEff <- resEff[!(resEff$resources < mean(resEff$resources) - sd(resEff$resources)*3 | resEff$resources > mean(resEff$resources)+ sd(resEff$resources)*3),]
  resEff <- merge(resEff, QA_df, by='scanid')
  resEff <- resEff[which(resEff[,2]>0),]
  if (i==999)
  {
    resEff$probability <- rep(i/1000,length(resEff[,1]))
  } else if (i>90 & i <100) {
    resEff$probability <- rep(i/100,length(resEff[,1]))
  } else {
    resEff$probability <- rep(i/10,length(resEff[,1]))
  }
  resEff <- as.data.frame(cbind(resEff$scanid,resEff$resources,resEff$pathComplexity,resEff$probability))
  colnames(resEff)<- c('scanid',"resources","pathComplexity","probability")
  assign(paste0("df",i),as.data.frame(resEff))
}

df1 <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df92,df94,df96,df98,df999)

df1$probability <- as.factor(df1$probability)

df2<- subset(df1,df1$probability==0.999)
df2 <- merge(df2, QA_df, by=c('scanid'))

gam1 <- gam(resources~ s(pathComplexity.x, k=4)+degree+density+s(ageAtScan1,k=4) +sex+s(ageAtScan1, by=sex, k=4) +meanMotion,fx=TRUE,method="REML",data=df2)

gam2 <- gam(resources~ pathComplexity.x+degree+density+s(ageAtScan1,k=4) +sex+s(ageAtScan1, by=sex, k=4) +meanMotion,fx=TRUE,method="REML",data=df2)

# Compare non-linear fit with linear fit: lower scores indicate better model
AIC(gam1)-AIC(gam2)
BIC(gam1)-BIC(gam2)

visreg(gam1,"pathComplexity.x",gg=T)+theme_classic(base_size=20)

ggsave('figures/fig6/lowFidelityRegime999.eps',device='eps',width=6.31,height=6.31)