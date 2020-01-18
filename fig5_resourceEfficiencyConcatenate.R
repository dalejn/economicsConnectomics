rm(list=ls())

################################
# Load workspace and libraries #
################################

load('data/0_finalData.RData')

########################################
# Concatenate resource efficiency send #
########################################

resEff_regional <- as.data.frame(cbind(rep(NA,360)))
column_index = 1
for (i in c(999,98,96,94,92,9,8,7,6,5,4,3,2,1))
{
  filepath <- paste0('data/resource_efficiency0',i,'_send.txt')
  
  resEff <- read.csv(filepath,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]
  
  resEffReg <- cbind(resEff[1],resEff[2:ncol(resEff)])
  is.na(resEffReg) <- do.call(cbind,lapply(resEffReg, is.infinite))
  resEff[resEff==0] <- NA
  
  column_new <- paste0('slope',seq(1:360))
  colnames(resEffReg) <- c('scanid',column_new)
  resEffReg <- merge(resEffReg, QA_df, by='scanid')
  
  for (j in 2:(length(column_new)+1)){
    resEff_regional[j-1,column_index] <- mean(resEffReg[[j]],na.rm=T)
    print(paste0(i,'_',j))
  }
  column_index = column_index+1
}

fileOutName <- paste0('data/resourceEfficiency_regional_send.csv')
write.table(resEff_regional,fileOutName,sep=',',row.names=F)

###########################################
# Concatenate resource efficiency receive #
###########################################

resEff_regional <- as.data.frame(cbind(rep(NA,360)))
column_index = 1
for (i in c(999,98,96,94,92,9,8,7,6,5,4,3,2,1))
{
  filepath <- paste0('data/resource_efficiency0',i,'_receive.txt')

  resEff <- read.csv(filepath,header=F,sep=' ')
  resEff[2:ncol(resEff)] <- 1/resEff[2:ncol(resEff)]

  resEffReg <- cbind(resEff[1],resEff[2:ncol(resEff)])
  is.na(resEffReg) <- do.call(cbind,lapply(resEffReg, is.infinite))
  resEff[resEff==0] <- NA

  column_new <- paste0('slope',seq(1:360))
  colnames(resEffReg) <- c('scanid',column_new)
  resEffReg <- merge(resEffReg, QA_df, by='scanid')

  for (j in 2:(length(column_new)+1)){
    resEff_regional[j-1,column_index] <- mean(resEffReg[[j]],na.rm=T)
    print(paste0(i,'_',j))
  }
  column_index = column_index+1
}

fileOutName <- paste0('data/resourceEfficiency_regional_receive.csv')
write.table(resEff_regional,fileOutName,sep=',',row.names=F)
