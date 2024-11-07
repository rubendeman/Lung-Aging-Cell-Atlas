setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

#read glmmTMB output
cell.types<-c('AT2B','AT2S')
names(cell.types)<-cell.types
#datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_AGE_GENES/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(datal,FUN = setNames,nm=c('Beta0','Beta1','p_val','p_val_adj'))
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})

noquote(datal[['AT2B']]$Gene[datal[['AT2B']]$Beta1<0&datal[['AT2B']]$p_val_adj<0.1])
