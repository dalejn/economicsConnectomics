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

# load global efficiency per individual data

globEff<- read.csv("/data/jux/BBL/projects/ASLnetwork/results/global_efficiency.mat",sep=" ",header=F)
colnames(globEff)<- c('scanid','globEff')

QA_df <- merge(QA_df,globEff,by=c("scanid"))
QA_df$ageInYears <- QA_df$ageAtScan1/12

##############################
# Global efficiency analyses #
##############################

# Assess global efficiency with development
globEffByAge<- gam(globEff ~ s(ageInYears,k=4) + degree + sex + s(ageInYears,k=4,by=sex)+ meanMotion,fx=TRUE,method="REML",data=QA_df)

visreg(globEffByAge,"ageInYears", xlab="Age (years)", ylab="Global Efficiency", gg=T)+theme_classic(base_size=20)+labs(y="Global efficiency", x="Age (years)")
ggsave('figures/fig2/globalEffAge.eps',device='eps',width=7.18,height=6.31)

summary(globEffByAge)

# Assess global CBF with development

# age effect
globEffByAge<- gam(globalCBF ~ s(ageInYears,k=4) + degree + sex +  meanMotion,fx=TRUE,method="REML",data=QA_df)
summary(globEffByAge)

# age-by-sex interaction (replicating Satterthwaite et al. 2014 PNAS)
globEffByAge<- gam(globalCBF ~ s(ageInYears,k=4) + degree + sex + s(ageInYears,k=4,by=sex)+  meanMotion,fx=TRUE,method="REML",data=QA_df)

visreg(globEffByAge,"ageInYears", xlab="Age (years)", ylab="Global CBF (ml/100g/min)", gg=T)+theme_classic(base_size=20)+labs(y="Global CBF (ml/100g/min)", x="Age (years)")
ggsave('figures/fig2/globalCbfAge.eps',device='eps',width=7.18,height=6.31)
summary(globEffByAge)

# Assess partial correlation of global efficiency and CBF, controlling for age

globEffByAge<- gam(globEff ~ s(ageAtScan1,k=4), method="REML",data=QA_df)
residGlobEff <- resid(globEffByAge)
globEffByAge<- gam(globalCBF ~ s(ageAtScan1,k=4), method="REML",data=QA_df)
residGlobCbf <- resid(globEffByAge)
cor.test(residGlobEff,residGlobCbf)

####################################
# Assess walk lengths individually #
####################################

# reload data
rm(list=ls())
load('0_finalData.RData')

l <- vector("list",15)
t <- vector("list",15)
degreesFreedom <- vector("list",15)

for (i in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
{
  filepath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/faExp/fa_exp',i,'_individual.txt')
  fa <- as.data.frame(read.table(filepath,header=F,sep=' '))
  column_new <- paste0('fa',i)
  colnames(fa) <- c('scanid',column_new)
  fa[[column_new]]<- scale(fa[[column_new]])
  QA_df<- merge(QA_df,fa,by=c('scanid'))
  
  placeholder <- QA_df[[column_new]]

  faCbf<- gam(globalCBF ~ as.vector(placeholder) + s(ageAtScan1,k=4) + degree + density + sex + meanMotion+ s(ageAtScan1,k=4,by=sex),fx=TRUE,method="REML",data=QA_df)
  l[[i]] <- summary(faCbf)$p.pv[2]
  t[[i]]<- summary(faCbf)$p.t[2]
  degreesFreedom[[i]] <- sum(faCbf$edf)
}
df1<- as.data.frame(do.call(rbind,l))
df2<- as.data.frame(do.call(rbind,t))
df3<- as.data.frame(do.call(rbind,degreesFreedom))

# Plot results

x_pathLength <- c(seq(1,15,1))
pvalues <- as.data.frame(c(df1[1]))
indPvals <- pvalues

colnames(pvalues)<- c('z-score')
y_zscore <- abs(qnorm(pvalues$'z-score', lower.tail=TRUE))

df1 <- data.frame(x_pathLength, y_zscore)
colnames(df1)<- c("Correlation of (FA^n) and CBF", "zscore")

library(grid)
p <-ggplot(df1,aes(df1[,1],df1[,2]))
p <- p + geom_col()
p <- p + expand_limits(x = 0, y = 0)
p <- p + xlim(0,16) +ylim(0,3)
p <- p +
  theme_minimal(base_size=30) + # start with a minimal theme and add what we need
  theme(text = element_text(color = "gray20",size=30),
        legend.position = c("top"), # position the legend in the upper left 
        legend.direction = "horizontal",
        legend.justification = 0.1, # anchor point for legend.position.
        legend.text = element_text(size = 17, color = "gray10"),
        axis.text = element_text(face = "italic", size= 25),
        axis.title.x = element_text(vjust = 0, size = 25), # move title away from axis
        axis.title.y = element_text(vjust = .5, size = 25), # move away for axis
        axis.ticks.y = element_blank(), # element_blank() is how we remove elements
        # axis.line = element_line(color = "gray40", size = 0.5),
        axis.line.y = element_blank(),
        panel.grid.major = element_line(color = "gray50", size = 0.5),
        panel.grid.major.x = element_blank()
  ) 
p <- p + labs(x=bquote('Walk length n in'~FA^n), y='z-scored relationship between \n global strength and global CBF')
p <- p + geom_hline(yintercept=abs(qnorm(0.05, lower.tail=TRUE)), color="royalblue", linetype="dashed",size=1)
p

ggsave('figures/fig2/interIndividualWalkLengthCbf.eps',device='eps',width=7.18,height=6.31)

########################################
# Assess walk lengths brain regionally #
########################################

# reload data
load('0_finalData.RData')

# load regional CBF
# verify regional CBF of Glasser parcels

colnames(QA_df[297:656])
df1 <- QA_df[297:656]
df1 <- df1[c(343,278,352,344,346,350,355,200,199,250,294,209,347,320,315,354,268,251,348,272,273,310,279,234,318,324,292,330,224,221,314,190,345,246,196,211,212,189,210,191,192,219,323,215,220,223,222,271,356,277,181,188,198,214,216,218,286,197,237,288,236,247,228,287,184,206,226,225,232,227,233,182,229,202,203,205,240,217,260,261,263,262,290,204,242,230,231,185,235,183,186,187,280,207,270,213,259,322,201,283,281,282,208,319,296,313,331,255,276,309,245,243,252,254,253,299,244,249,316,258,317,293,325,291,241,304,306,326,327,329,338,332,334,335,337,336,307,303,340,341,342,248,301,267,266,265,298,295,297,300,302,353,357,359,305,351,257,349,274,358,195,194,360,193,321,311,312,264,256,284,289,339,275,269,239,328,333,308,238,285,163,98,172,164,166,170,175,20,19,70,114,29,167,140,135,174,88,71,168,92,93,130,99,54,138,144,112,150,44,41,134,10,165,66,16,31,32,9,30,11,12,39,143,35,40,43,42,91,176,97,1,8,18,34,36,38,106,17,57,108,56,67,48,107,4,26,46,45,52,47,53,2,49,22,23,25,60,37,80,81,83,82,110,24,62,50,51,5,55,3,6,7,100,27,90,33,79,142,21,103,101,102,28,139,116,133,151,75,96,129,65,63,72,74,73,119,64,69,136,78,137,113,145,111,61,124,126,146,147,149,158,152,154,155,157,156,127,123,160,161,162,68,121,87,86,85,118,115,117,120,122,173,177,179,125,171,77,169,94,178,15,14,180,13,141,131,132,84,76,104,109,159,95,89,59,148,153,128,58,105)]
df1 <- as.data.frame(cbind(QA_df[1],df1))

roiCbf <- cbind(colMeans(df1[2:361],na.rm=T))

degreeRoi <- unlist(read.csv('/data/jux/BBL/projects/ASLnetwork/results/degreeRoi.txt',header=F,sep=' '))

r <- vector("list",15)
p <- vector("list",15)

for (i in c(seq(1,15,1)))
{
  filepath <- paste0('/data/jux/BBL/projects/ASLnetwork/results/faExp/fa_exp',i,'_regional.txt')
  fa <- read.table(filepath,header=F,sep=' ')
  column_new <- paste0('fa',seq(1:360))
  colnames(fa) <- c('scanid',column_new)
  
  for (j in column_new){
      fa[[j]]<- scale(fa[[j]])
  }

  fa <- merge(fa, QA_df, by='scanid')
  
  for (j in column_new){
    lm1 <- gam(as.vector(fa[[j]]) ~ s(ageAtScan1,by=sex,k=4) + s(ageAtScan1,k=4) + sex + degree + density + meanMotion,fx=TRUE,method="REML",data=fa)
    fa[j] <- resid(lm1)
    
    print(paste0(i,'_',j))
  }
  
  roiFa <- cbind(colMeans(fa[,2:361],na.rm=T))
  roiFa_resid <- roiFa
  r[[i]]<-cor.test(roiCbf,roiFa_resid,method="spearman")$estimate
  p[[i]]<-cor.test(roiCbf,roiFa_resid,method="spearman")$p.value
}

df1<- as.data.frame(do.call(rbind,r))
df2 <- as.data.frame(do.call(rbind,p))
regionalPvals <- df2

x_pathLength <- c(seq(1,15,1))
pvalues <- as.data.frame(c(df2[1]))
colnames(pvalues)<- c('z-score')
y_zscore <- pvalues
y_zscore[pvalues$`z-score` < 0.5,] <- abs(qnorm(pvalues$'z-score'[pvalues$`z-score` <0.5], lower.tail=TRUE))
y_zscore[pvalues$`z-score` > 0.5,] <- 0

df1 <- data.frame(x_pathLength, y_zscore)

library(grid)
p <-ggplot(df1,aes(df1$x_pathLength,df1$z.score))
p <- p + geom_col()
p <- p + expand_limits(x = 0, y = 0)
p <- p + xlim(0,16) +ylim(0,3)
p <- p +
  theme_minimal(base_size=25) + # start with a minimal theme and add what we need
  theme(text = element_text(color = "gray20"),
        legend.position = c("top"), # position the legend in the upper left 
        legend.direction = "horizontal",
        legend.justification = 0.1, # anchor point for legend.position.
        legend.text = element_text(size = 17, color = "gray10"),
        axis.text = element_text(face = "italic", size= 25),
        axis.title.x = element_text(vjust = 0, size = 25), # move title away from axis
        axis.title.y = element_text(vjust = .5, size = 25), # move away for axis
        axis.ticks.y = element_blank(), # element_blank() is how we remove elements
        # axis.line = element_line(color = "gray40", size = 0.5),
        axis.line.y = element_blank(),
        panel.grid.major = element_line(color = "gray50", size = 0.5),
        panel.grid.major.x = element_blank()
  )

p <- p + geom_hline(yintercept=abs(qnorm(0.05, lower.tail=TRUE)), color="royalblue", linetype="dashed",size=1)
p <- p + labs(x=bquote('Walk length'~FA^n), y='z-scored relationship between \n regional strength and regional CBF')
p

ggsave('figures/fig2/regionalWalkLengthCbf.eps',device='eps',width=7.18,height=6.31)

##################
# FDR correction #
##################

allPvals <- as.vector(unlist(c(indPvals,pvalues)))
p.adjust(allPvals, method='fdr')