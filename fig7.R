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

####################################

# load rich club assignments
richclub <- read.csv('/data/jux/BBL/projects/ASLnetwork/results/richclub.txt')
df1<- cbind(richclub)

# For right stochastic matrices:
# load compression efficiency send as row mean of resource efficiency over i
# load compression efficiency receive as column mean of resource efficiency over j

#unbiasedSlopes_send <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/tSpectrum/distortion_regional.txt', sep=' ',header=F))
#unbiasedSlopes_receive <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/minus3distortion_regional.txt', sep=' ',header=F))

unbiasedSlopes_send <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_send.txt', sep=' ',header=F))
unbiasedSlopes_receive <- as.data.frame(read.csv('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_receive.txt', sep=' ',header=F))

df1$slope_send <- unbiasedSlopes_send[,2]
df1$slope_receive <- unbiasedSlopes_receive[,2]

scaling <- as.vector(read.csv('/data/jux/BBL/projects/ASLnetwork/results/regionalData/allometricScaling_regional.txt', sep=' ',header=F))
df1$scaling <- scaling$V1

myelin <- as.vector(read.csv('/data/jux/BBL/projects/ASLnetwork/results/regionalData/myelin_regional.txt', sep=' ',header=F))
df1$myelin <- myelin$V1

colnames(df1) <- c('richclub','slope_send', 'slope_receive', 'scaling', 'myelin')

####################################
# Uncorrected spatial correlations #
####################################

myelinP <- cor.test(df1$myelin,df1$slope_receive, method='spearman')$p.value
myelinP2 <- cor.test(df1$myelin,df1$slope_send, method='spearman')$p.value

scalingP <- cor.test(df1$scaling,df1$slope_receive, method='spearman')$p.value
scalingP2 <- cor.test(df1$scaling,df1$slope_send, method='spearman')$p.value

# uncorrected p-values
evo_pvals <- c(scalingP2,myelinP2, scalingP, myelinP)
p.adjust(evo_pvals, method='holm')

###################################
# Spatially permuted correlations #
###################################

# see Matlab code fig7_atlasCorrelations.m

evo_spin_pvals <- c(0.0512, 0.0097, 0.0069, 0.0226)
p.adjust(evo_spin_pvals, method='holm')

#####################
# Plot correlations #
#####################

ggplot(df1,aes(slope_send, myelin))+geom_point()+geom_smooth(method='lm')+theme_classic(base_size=20)
ggsave("figures/fig7/myelin-distortion_send.eps",device='eps',width=7.18,height=6.31)

ggplot(df1, aes(slope_receive, scaling))+geom_point()+geom_smooth(method='lm')+theme_classic(base_size=20)
ggsave("figures/fig7/scaling-distortion_receive.eps",device='eps',width=7.18,height=6.31)

ggplot(df1, aes(slope_send, scaling))+geom_point()+geom_smooth(method='lm')+theme_classic(base_size=20)
ggsave("figures/fig7/scaling-distortion_send.eps",device='eps',width=7.18,height=6.31)