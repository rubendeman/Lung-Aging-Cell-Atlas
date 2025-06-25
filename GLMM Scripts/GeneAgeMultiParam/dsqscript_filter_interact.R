.libPaths(c("/home/rd796/project/R/4.2", .libPaths()))
library(Seurat)
library(tidyverse)
library(dplyr)

####### set working directory
my.workingDir <- "/home/rd796/palmer_scratch"

setwd(my.workingDir)

load('/home/rd796/project/ageproj/agingseurat.RData')
DefaultAssay(immune.combined)<-'RNA'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'

#comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
#cell.types<-unlist(comp)
cell.types<-c('AT1','AT2','gCap')

genes.test <- list() #keep genes expressed in greater than 50% of subjects
for (i in cell.types){
tmp=subset(immune.combined,subset=predicted.id==i);
countdata=tmp@assays$RNA@counts
genekeep<-rownames(tmp)[rowSums(countdata>0)>ncol(countdata)*0.05]
sids=unique(tmp$orig.ident)
sidtbl<-lapply(sids,function(x){countdata[,tmp$orig.ident==x]})
sidtbl<-sidtbl[!sapply(sapply(sidtbl,FUN=ncol),FUN=is.null)]
sidtbl2<-lapply(sidtbl,function(y){rownames(tmp)[rowSums(y>0)>0]}) #ncol(y)*0.01
sumtbl<-as.data.frame(table(unlist(sidtbl2)))
sumtbl$Var1<-as.character(sumtbl$Var1)
genes.test[[i]]<-sumtbl$Var1[sumtbl$Freq>length(sidtbl2)*0.5&sumtbl$Var1 %in% genekeep]
}

#### we need to split the genes up across jobs
master.gene.groups <- lapply(genes.test,function(z){split(z, ceiling(seq_along(z)/100))})

#### define where all the scripts are going to be written to
script.output.dir <- file.path(my.workingDir,"glmm.scripts")
### make this folder if it doesn't already exist
if(dir.exists(script.output.dir)==FALSE){
  dir.create(script.output.dir)
}

jobsub.filepath <- file.path(script.output.dir, "multiparam_joblist.txt")

rscript.function.filepath <- "/gpfs/gibbs/pi/kaminski/public/Backup/Ruben/GLMMAging/GeneAgeMultiParam/runmodel_interact.R"

### what cell type to test
cellType=cell.types

### create a new jobsub file before starting
cat("", sep="", file=jobsub.filepath, append=FALSE)

cmd.out <- NULL
for(k in c('Sex','Dataset','Smoke')){
for(j in 1:length(cellType)){
gene.groups=master.gene.groups[[j]]
for(i in 1:length(gene.groups)){
  temp.param <- k
  temp.cellType <- cellType[j]
  temp.genes <- gene.groups[[i]] 
  cmd.out <- paste("module load R/4.2.0-foss-2020b; ", sep="")
  cmd.out <- paste(cmd.out, "Rscript ", rscript.function.filepath, " ", paste(c(temp.param, temp.cellType, temp.genes), collapse=" "), sep="")
  cat(cmd.out, "\n", sep="", file=jobsub.filepath, append=TRUE)
  }
}
}