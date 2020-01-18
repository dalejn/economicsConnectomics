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

# load rich club assignments
richclub <- read.csv('/data/richclub.txt')
df1<- cbind(richclub)

############################
# Rich club CBF Comparison #
############################

# fix the order of regions which got stored in the dataframe
# as alphabetical. Change back to Glasser ordering

# load CBF data (order will need to be unscrambled from alphabetical sorting in data frame)

# print scrambled order of brain regions
colnames(QA_df[297:656])
roiCbf <- QA_df[297:656]

# unscramble brain regions using the Glasser indices
roiCbf <- roiCbf[c(343,278,352,344,346,350,355,200,199,250,294,209,347,320,315,354,268,251,348,272,273,310,279,234,318,324,292,330,224,221,314,190,345,246,196,211,212,189,210,191,192,219,323,215,220,223,222,271,356,277,181,188,198,214,216,218,286,197,237,288,236,247,228,287,184,206,226,225,232,227,233,182,229,202,203,205,240,217,260,261,263,262,290,204,242,230,231,185,235,183,186,187,280,207,270,213,259,322,201,283,281,282,208,319,296,313,331,255,276,309,245,243,252,254,253,299,244,249,316,258,317,293,325,291,241,304,306,326,327,329,338,332,334,335,337,336,307,303,340,341,342,248,301,267,266,265,298,295,297,300,302,353,357,359,305,351,257,349,274,358,195,194,360,193,321,311,312,264,256,284,289,339,275,269,239,328,333,308,238,285,163,98,172,164,166,170,175,20,19,70,114,29,167,140,135,174,88,71,168,92,93,130,99,54,138,144,112,150,44,41,134,10,165,66,16,31,32,9,30,11,12,39,143,35,40,43,42,91,176,97,1,8,18,34,36,38,106,17,57,108,56,67,48,107,4,26,46,45,52,47,53,2,49,22,23,25,60,37,80,81,83,82,110,24,62,50,51,5,55,3,6,7,100,27,90,33,79,142,21,103,101,102,28,139,116,133,151,75,96,129,65,63,72,74,73,119,64,69,136,78,137,113,145,111,61,124,126,146,147,149,158,152,154,155,157,156,127,123,160,161,162,68,121,87,86,85,118,115,117,120,122,173,177,179,125,171,77,169,94,178,15,14,180,13,141,131,132,84,76,104,109,159,95,89,59,148,153,128,58,105)]
roiCbf <- as.data.frame(cbind(QA_df[1],roiCbf))

#df1$roiCbf <- as.numeric(cbind(colMeans(roiCbf[2:361],na.rm=T)))

column_new <- paste0('cbf',seq(1:360))
colnames(roiCbf) <- c('scanid',column_new)
roiCbf <- merge(roiCbf, QA_df, by='scanid')
roiCbf_sensitivity <- merge(roiCbf, QA_df, by='scanid')
for (j in column_new){
  lm1 <- gam(as.vector(roiCbf[[j]]) ~ s(ageAtScan1,by=sex,k=4) + s(ageAtScan1,k=4)+ sex + meanMotion,fx=TRUE,method="REML", na.action = na.exclude, data=roiCbf)
  lm2 <- gam(as.vector(roiCbf[[j]]) ~ s(ageAtScan1,by=sex,k=4) + s(ageAtScan1,k=4)+ sex,fx=TRUE,method="REML", na.action = na.exclude, data=roiCbf)
  roiCbf[j] <- resid(lm1)
  roiCbf_sensitivity[j] <- resid(lm2)
  print(paste0(colnames(roiCbf[j])))
}

roiCbf <- cbind(colMeans(roiCbf[,2:361],na.rm=T))
roiCbf_sensitivity <- cbind(colMeans(roiCbf_sensitivity[,2:361],na.rm=T))

# no covariates regional
# df1 <- as.data.frame(cbind(as.numeric(cbind(colMeans(roiCbf[2:361],na.rm=T))), richclub))

# with covariates regional
df1 <- as.data.frame(cbind(as.numeric(roiCbf), richclub))
colnames(df1) <- c('cbf','richclub')
df1$richclub <- as.factor(df1$richclub)

df2 <- as.data.frame(cbind(as.numeric(roiCbf_sensitivity), richclub))
colnames(df2) <- c('cbf','richclub')
df2$richclub <- as.factor(df2$richclub)

#rich club CBF versus non-richclub CBF
wilcox.test(df1$cbf[which(df1$richclub==0)], df1$cbf[which(df1$richclub==1)])

wilcox.test(df2$cbf[which(df2$richclub==0)], df2$cbf[which(df2$richclub==1)])

# difference between rich club and non rich club CBF
ggplot(data=df1,aes(x=df1$richclub,y=df1$cbf))+geom_boxplot(outlier.size=0.05)+theme_classic(base_size=20)
ggsave('figures/supplement/richclubCBF.eps',device='eps',width=7.18,height=6.31)

ggplot(data=df1,aes(x=df1$richclub,y=df1$cbf))+geom_boxplot(outlier.size=0.05)+theme_classic(base_size=20)
ggsave('figures/supplement/richclubCBF_sensitivity.eps',device='eps',width=7.18,height=6.31)
