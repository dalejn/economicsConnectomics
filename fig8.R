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

# load rich club assignments
richclub <- read.csv('/data/jux/BBL/projects/ASLnetwork/results/richclub.txt')
df1<- cbind(richclub)

# For right stochastic matrices:
# load compression efficiency send as row mean of resource efficiency over i
# load compression efficiency receive as column mean of resource efficiency over j

unbiasedSlopes_send <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_send.txt', sep=' ',header=F))
unbiasedSlopes_receive <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_receive.txt', sep=' ',header=F))

unbiasedSlopes_send <- as.data.frame(unbiasedSlopes_send[,2])
unbiasedSlopes_receive <- as.data.frame(unbiasedSlopes_receive[,2])

####################################
# Rich club compression efficiency #
####################################

unbiasedSlopes_receive$direction <- rep('receive', 360)
colnames(unbiasedSlopes_receive) <- c('slope','direction')
unbiasedSlopes_send$direction <- rep('send', 360)
colnames(unbiasedSlopes_send) <- c('slope','direction')
combined_slopes <- rbind(unbiasedSlopes_receive, unbiasedSlopes_send)

df1 <- as.data.frame(cbind(combined_slopes, richclub))
colnames(df1) <- c('slopes','direction', 'richclub')
df1$richclub <- as.factor(df1$richclub)
df1$direction <- as.factor(df1$direction)

setEPS()
postscript("figures/fig8/slopeByRichClub_sendReceive.eps",width=14.18,height=6.31)
ggplot(data=df1,aes(df1$direction,df1$slopes,fill=df1$richclub))+geom_boxplot(outlier.size=2)+theme_classic(base_size=20)
dev.off()

wilcox.test(df1$slope[df1$direction=="send" & df1$richclub==0], df1$slope[df1$direction=="send" & df1$richclub==1])

wilcox.test(df1$slope[df1$direction=="receive" & df1$richclub==0], df1$slope[df1$direction=="receive" & df1$richclub==1])

#######################################################################
# Correlation between compression efficiency and cognitive efficiency #
#######################################################################

# load individual compression efficiency measures (averaged across rich-club or non-richclub)

unbiasedSlopes <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/slope_unbiased_richclub0.txt', sep=' ',header=F))
colnames(unbiasedSlopes) <- c('scanid','unbiasedSlopes_richclub0')
QA_df <- merge(QA_df, unbiasedSlopes, by=c('scanid'))

unbiasedSlopes <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/slope_unbiased_richclub1.txt', sep=' ',header=F))
colnames(unbiasedSlopes) <- c('scanid','unbiasedSlopes_richclub1')
QA_df <- merge(QA_df, unbiasedSlopes, by=c('scanid'))


# load individual cognitive efficiency measures

QA_df$F1_Complex_Reasoning_Efficiency <- as.numeric(as.character(QA_df$F1_Complex_Reasoning_Efficiency))
QA_df$F2_Memory_Efficiency <- as.numeric(as.character(QA_df$F2_Memory.Efficiency))
QA_df$F3_Executive_Efficiency <- as.numeric(as.character(QA_df$F3_Executive_Efficiency))
QA_df$F4_Social_Cognition_Efficiency <- as.numeric(as.character(QA_df$F4_Social_Cognition_Efficiency))

# Exclude outliers outside 3 SDs

QA_minusOutliers <- QA_df[which(QA_df$unbiasedSlopes_richclub1 > mean(QA_df$unbiasedSlopes_richclub1)-sd(QA_df$unbiasedSlopes_richclub1)*3),]
QA_minusOutliers <- QA_minusOutliers[which(QA_minusOutliers$unbiasedSlopes_richclub1 < mean(QA_minusOutliers$unbiasedSlopes_richclub1)+sd(QA_minusOutliers$unbiasedSlopes_richclub1)*3),]
QA_minusOutliers <- QA_minusOutliers[which(QA_minusOutliers$unbiasedSlopes_richclub0 > mean(QA_minusOutliers$unbiasedSlopes_richclub0)-sd(QA_minusOutliers$unbiasedSlopes_richclub0)*3),]
QA_minusOutliers <- QA_minusOutliers[which(QA_minusOutliers$unbiasedSlopes_richclub0 < mean(QA_minusOutliers$unbiasedSlopes_richclub0)+sd(QA_minusOutliers$unbiasedSlopes_richclub0)*3),]

# load individual global efficiency measures

globEff<- read.csv("/data/jux/BBL/projects/ASLnetwork/results/global_efficiency.mat",sep=" ",header=F)
colnames(globEff)<- c('scanid','globEff')
QA_df <- merge(QA_df,globEff,by=c("scanid"))

# Compression efficiency vs cognitive efficiency for non-richclub

    globEffCbf<- gam(F1_Complex_Reasoning_Efficiency ~ unbiasedSlopes_richclub0  + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    complexReasoningEfficiencyP <- summary(globEffCbf)$p.pv[2]
    complexReasoningEfficiencyT <- summary(globEffCbf)$p.t[2]    
    complexReasoningEfficiencyDf <- sum(globEffCbf$edf)

    globEffCbf<- gam(F2_Memory_Efficiency ~ unbiasedSlopes_richclub0 + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    memoryEfficiencyP <- summary(globEffCbf)$p.pv[2]
    memoryEfficiencyT <- summary(globEffCbf)$p.t[2]
    memoryEfficiencyDf <- sum(globEffCbf$edf)
     
    globEffCbf<- gam(F3_Executive_Efficiency ~ unbiasedSlopes_richclub0 + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df)    
    execEfficiencyP <- summary(globEffCbf)$p.pv[2]
    execEfficiencyT <- summary(globEffCbf)$p.t[2]    
    execEfficiencyDf <- sum(globEffCbf$edf)

    globEffCbf<- gam(F4_Social_Cognition_Efficiency ~ unbiasedSlopes_richclub0 + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    socialCognitionEfficiencyP <- summary(globEffCbf)$p.pv[2]
    socialCognitionEfficiencyT <- summary(globEffCbf)$p.t[2]    
    socialCognitionEfficiencyDf <- sum(globEffCbf$edf)

    NCB_pvals <- c(complexReasoningEfficiencyP,memoryEfficiencyP,execEfficiencyP,socialCognitionEfficiencyP)
    NCB_tvals <- c(complexReasoningEfficiencyT,memoryEfficiencyT,execEfficiencyT,socialCognitionEfficiencyT)
    NCB_df <- c(complexReasoningEfficiencyDf,memoryEfficiencyDf,execEfficiencyDf,socialCognitionEfficiencyDf)

# Compression efficiency vs cognitive efficiency for rich-club
    
    globEffCbf<- gam(F1_Complex_Reasoning_Efficiency ~ unbiasedSlopes_richclub1 + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    complexReasoningEfficiencyP <- summary(globEffCbf)$p.pv[2]
    complexReasoningEfficiencyT <- summary(globEffCbf)$p.t[2]
    complexReasoningEfficiencyDf <- sum(globEffCbf$edf)
    

    globEffCbf<- gam(F2_Memory_Efficiency ~ unbiasedSlopes_richclub1+ s(ageAtScan1,k=4)  + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    memoryEfficiencyP <- summary(globEffCbf)$p.pv[2]
    memoryEfficiencyT <- summary(globEffCbf)$p.t[2]
    memoryEfficiencyDf <- sum(globEffCbf$edf)    

    globEffCbf<- gam(F3_Executive_Efficiency ~ unbiasedSlopes_richclub1 +s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df)
    execEfficiencyP <- summary(globEffCbf)$p.pv[2]
    execEfficiencyT <- summary(globEffCbf)$p.t[2]
    execEfficiencyDf <- sum(globEffCbf$edf)    

    globEffCbf<- gam(F4_Social_Cognition_Efficiency ~ unbiasedSlopes_richclub1 + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_minusOutliers)
    socialCognitionEfficiencyP <- summary(globEffCbf)$p.pv[2]
    socialCognitionEfficiencyT <- summary(globEffCbf)$p.t[2]
    socialCognitionEfficiencyDf <- sum(globEffCbf$edf)

    NCB_pvals2 <- c(complexReasoningEfficiencyP,memoryEfficiencyP,execEfficiencyP,socialCognitionEfficiencyP)
    NCB_tvals <- c(complexReasoningEfficiencyT,memoryEfficiencyT,execEfficiencyT,socialCognitionEfficiencyT)
    NCB_df <- c(complexReasoningEfficiencyDf,memoryEfficiencyDf,execEfficiencyDf,socialCognitionEfficiencyDf)

pvals <- c(NCB_pvals,NCB_pvals2)
p.adjust(pvals,method='holm')

######################################################################
# Plot histograms of compression efficiency and cognitive efficiency #
######################################################################

## Calculate Z-scores 

NCB_Zscores <- c(qnorm(NCB_pvals, lower.tail=FALSE),qnorm(NCB_pvals2, lower.tail=FALSE))
NCB_Zscores <- as.data.frame(NCB_Zscores)
NCB_Zscores$System_idx <- c(1:4,1:4)

NCB_Zscores$Effect <- base::rank(-NCB_Zscores$NCB_Zscores)
NCB_Zscores$richclub <- c(0,0,0,0,1,1,1,1)

p <- ggplot(data=NCB_Zscores, aes(System_idx, NCB_Zscores,fill=as.factor(richclub))) + geom_bar(col="black",stat="identity",position=position_dodge())
p<- p + scale_x_continuous(breaks=1:4, labels=c("Complex Reasoning","Memory","Executive","Social Cognition"))+ theme(axis.text = element_text(size= 20), axis.text.x=element_text(angle = 60, hjust = 1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border= element_rect(colour="black",fill=FALSE, size=1))

ggsave('figures/fig8/NCB_zscore_byChannelCapacity.eps',device='eps',width=7.18,height=6.31)

##################################################################
# Scatterplot of compression efficiency and cognitive efficiency #
##################################################################
QA_df$unbiasedSlopes <- NULL
# Load individual measure of compression efficiency

unbiasedSlopes <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/unbiasedSlopes.txt', sep=' ',header=F))
unbiasedSlopes[,2] <- scale(unbiasedSlopes[,2])
colnames(unbiasedSlopes) <- c('scanid','unbiasedSlopes')
QA_df <- merge(QA_df, unbiasedSlopes, by=c('scanid'))

# checking collinearity of compression efficiency and global efficiency
globEff_resid <- resid(gam(globEff ~ s(ageAtScan1,k=4)+ sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df))
compEff_resid <- resid(gam(unbiasedSlopes ~ s(ageAtScan1,k=4)+ sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df))
cor.test(globEff_resid, compEff_resid, method="pearson")

# Controlling for global efficiecy

globEffCbf<- gam(F1_Complex_Reasoning_Efficiency ~ as.vector(unbiasedSlopes) + globEff + s(ageAtScan1,k=4) + sex+ s(ageAtScan1,by=sex,k=4) + degree+density + dti64MeanRelRMS,fx=TRUE,method="REML",data=QA_df)
summary(globEffCbf)

visreg(globEffCbf,"unbiasedSlopes",xlab="unbiased slopes",ylab="Complex Reasoning Efficiency",gg=T)+theme_classic(base_size=20)
ggsave('figures/fig8/complexReasoningEfficiencyByChannelCapacity.eps',device='eps',width=6.31,height=6.31)