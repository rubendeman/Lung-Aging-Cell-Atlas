setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(dplyr)
library(Seurat)

#LOAD INTEGRATED
load('imm12_10.RData')

DefaultAssay(immune.combined)<-'RNA'
immune.combined<-NormalizeData(immune.combined)

comp<-list(Epithelial=c('AT1','AT2','AT2B','AT2S','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage','DC'),Lymphoid=c('B','T','Mast','NK'))
cell.types<-unlist(comp)

immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)

comp_data<-immune.combined$predicted.id
for(i in names(comp)){comp_data[comp_data %in% comp[[i]]]<-i}
immune.combined<-AddMetaData(immune.combined,comp_data,col.name='comp')

#ADD METDATA
meta<-read.csv('demo2_17.csv')
sids<-immune.combined$orig.ident
inds<-match(sids,meta$Sample.ID)
meta=meta %>% mutate(agebin=cut(Age,breaks=c(0,60.5,100),labels=c('Young','Aged')))
meta=meta %>% mutate(agetri=cut(Age,breaks=c(0,40.5,60.5,100),labels=c('Young','Middle-Aged','Aged')))
meta2=meta[inds,]

immune.combined$age<-meta2$Age
immune.combined$sex<-meta2$Sex
immune.combined$smoke<-meta2$Smoking
immune.combined$agebin<-meta2$agebin
immune.combined$agetri<-meta2$agetri
immune.combined$Manuscript_Identity<-meta2$Dataset

save(immune.combined,file='agingseurat.RData')