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

# For right stochastic matrices:
# load compression efficiency send as row mean of resource efficiency over i
# load compression efficiency receive as column mean of resource efficiency over j

unbiasedSlopes_send <- as.data.frame(read.csv('/data/compressionEfficiency_send.txt', sep=' ',header=F))
unbiasedSlopes_receive <- as.data.frame(read.csv('/data/compressionEfficiency_receive.txt', sep=' ',header=F))

unbiasedSlopes_send <- as.data.frame(unbiasedSlopes_send[,2])
unbiasedSlopes_receive <- as.data.frame(unbiasedSlopes_receive[,2])

#####################################
# Subnetwork compression efficiency #
#####################################

# load data frame
df1 <- as.data.frame(cbind(unbiasedSlopes_send, unbiasedSlopes_receive, richclub))
colnames(df1) <- c('send','receive', 'richclub')
df1$richclub <- as.factor(df1$richclub)
subnetworks <- read.csv('/data/glasserInYeo.txt',header=F)
df1$subnetwork<- as.factor(subnetworks[,1])

# densities
ggplot(df1,aes(send,color=subnetwork))+geom_freqpoly(binwidth=0.05)+theme_classic(base_size=20)
ggsave("figures/supplement/send_distortionBySubnetwork_freqpoly.eps",device='eps',width=7.18,height=6.31)

ggplot(df1,aes(receive,color=subnetwork))+geom_freqpoly(binwidth=0.25)+theme_classic(base_size=20)
ggsave("figures/supplement/receive_distortionBySubnetwork_freqpoly.eps",device='eps',width=7.18,height=6.31)

# conditional probs
ggplot(data=df1,aes(send,stat(count),fill=subnetwork))+geom_density(position="fill")+theme_classic(base_size=20)
ggsave("figures/supplement/send_distortionBySubnetwork_proportion.eps",device='eps',width=7.18,height=6.31)

ggplot(data=df1,aes(receive,stat(count),fill=subnetwork))+geom_density(position="fill")+theme_classic(base_size=20)
ggsave("figures/supplement/receive_distortionBySubnetwork_proportion.eps",device='eps',width=7.18,height=6.31)

# histogram of means
subnetworkDistortion <- c(mean(df1$send[which(df1$subnetwork==1)]),mean(df1$send[which(df1$subnetwork==2)]),mean(df1$send[which(df1$subnetwork==3)]),mean(df1$send[which(df1$subnetwork==4)]),mean(df1$send[which(df1$subnetwork==5)]),mean(df1$send[which(df1$subnetwork==6)]),mean(df1$send[which(df1$subnetwork==7)]))
df1 <- as.data.frame(cbind(c(1:7), subnetworkDistortion))
df1$V1 <- as.factor(df1$V1)
ggplot(df1, aes(V1, subnetworkDistortion,fill=V1))+geom_bar(col="black",stat="identity")+theme_classic(base_size=20)
ggsave("figures/supplement/send_distortionBySubnetwork_histogram.eps",device='eps',width=7.18,height=6.31)

df1 <- as.data.frame(cbind(unbiasedSlopes_send, unbiasedSlopes_receive, richclub))
colnames(df1) <- c('send','receive', 'richclub')
df1$richclub <- as.factor(df1$richclub)
subnetworks <- read.csv('/data/glasserInYeo.txt',header=F)
df1$subnetwork<- as.factor(subnetworks[,1])

subnetworkDistortion <- c(mean(df1$receive[which(df1$subnetwork==1)]),mean(df1$receive[which(df1$subnetwork==2)]),mean(df1$receive[which(df1$subnetwork==3)]),mean(df1$receive[which(df1$subnetwork==4)]),mean(df1$receive[which(df1$subnetwork==5)]),mean(df1$receive[which(df1$subnetwork==6)]),mean(df1$receive[which(df1$subnetwork==7)]))
df1 <- as.data.frame(cbind(c(1:7), subnetworkDistortion))
df1$V1 <- as.factor(df1$V1)
ggplot(df1, aes(V1, subnetworkDistortion,fill=V1))+geom_bar(col="black",stat="identity")+theme_classic(base_size=20)
ggsave("figures/supplement/receive_distortionBySubnetwork_histogram.eps",device='eps',width=7.18,height=6.31)