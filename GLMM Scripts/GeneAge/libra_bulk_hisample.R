#.libPaths(c("/home/rd796/project/R/4.2", .libPaths()))
.libPaths(c("/home/rd796/project/R/4.2","/home/rd796/R/4.2", .libPaths()))
library(Seurat)
library(boot)
library(Libra)
library(tidyverse)
library(future)

plan("multicore", workers=2)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed=TRUE)

####### set working directory
my.workingDir <- "/home/rd796/palmer_scratch"

setwd(my.workingDir)

####### load the seurat object
load('/home/rd796/project/ageproj/imm12_10.RData')
DefaultAssay(immune.combined)<-'RNA'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'

comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)

immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)

#Hisamples only
immune.combined$temp<-paste0(immune.combined$predicted.id,"_",immune.combined$orig.ident)
tbl<-as.data.frame(table(immune.combined$temp))
keep<-tbl$Var1[tbl$Freq>20]
immune.combined<-subset(immune.combined,subset=temp %in% keep)

########### defining the parameters for yunqing's code
library(Libra)
meta<-immune.combined@meta.data
res.t <- run_de(immune.combined, meta=meta, cell_type_col="predicted.id",label_col="agebin",replicate_col="orig.ident")

### define where the output is being written to 
save(res.t,file='/home/rd796/palmer_scratch/libra_out_hisample.RData')