.libPaths(c("/home/rd796/project/R/4.2", .libPaths()))
library(Seurat)
library(boot)
library(glmmTMB)
library(tidyverse)
library(future)

plan("multicore", workers=2)
options(future.globals.maxSize = 5000 * 1024^2)
options(future.seed=TRUE)

####### set working directory
my.workingDir <- "/home/rd796/palmer_scratch/multi"

setwd(my.workingDir)

#### read in the arguments for the Rscript
inputArguments <- commandArgs(trailingOnly=TRUE)
paramTest <- inputArguments[1]
cellTypeTest <- inputArguments[2]
genes.test<- inputArguments[3:length(inputArguments)]

####### load the seurat object
l=load('/home/rd796/project/ageproj/agingseurat.RData')

DefaultAssay(immune.combined)<-'RNA'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='Alv. Fibroblast']='Alv_Fibroblast'
immune.combined$predicted.id[immune.combined$predicted.id=='Alv. Macrophage']='Alv_Macrophage'
immune.combined$predicted.id[immune.combined$predicted.id=='Adventitial Fibroblast']='Adv_Fibroblast'

comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)
immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)
immune.combined<-subset(immune.combined,subset=age>19,invert=F)
if (paramTest %in% immune.combined$sex){
    immune.combined<-subset(immune.combined,subset=sex == paramTest)
}else if (paramTest %in% c(1,2)){
    immune.combined<-subset(immune.combined,subset=Manuscript_Identity == paste("Dataset",paramTest))
}else if (paramTest %in% immune.combined$smoke){
immune.combined<-subset(immune.combined,subset=smoke == paramTest)
}else if(paramTest == 'Middle-Aged'){
    #immune.combined<-subset(immune.combined,subset=age<61)
    immune.combined$agebin[immune.combined$age>40]<-'Aged'
}else if (paramTest == 'Aged'){
    #immune.combined<-subset(immune.combined,subset=age>60|age<41)
}

soup.subset <- subset(immune.combined, subset=predicted.id==cellTypeTest)
### take all cells if below 1000
#nCells <- ncol(soup.subset)
#if(nCells > 5000){
#  set.seed(7)
#  cells.keep <- sample(colnames(soup.subset), size=5000, replace=FALSE)
#  soup.subset <- soup.subset[,cells.keep]
#}
gc()

########### defining the parameters for yunqing's code
gene <- genes.test
subject <- soup.subset@meta.data$orig.ident
sf <- soup.subset@meta.data$nCount_RNA
#age <- soup.subset@meta.data$age
age <- as.factor(soup.subset@meta.data$agebin)

### subset from the raw count matrix only the genes we are testing
count_mat <- soup.subset@assays$RNA@counts[gene, ]

### print to screen/sh.o what this job is gonna do
cat("Running GLMM\n", sep="")

############################################################### Ruben's code
res.t <- apply(count_mat, 1, function(gene){
  tmp.df <- data.frame(y=gene, norm_sf=sf, age=as.numeric(age), sub=as.numeric(factor(subject)))
  #f2 <- try(glmer(y ~ age + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = poisson))
  f2 <- try(glmmTMB(y ~ age + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = nbinom2)) #zi = ~ age;
    res2 <-try(c(
    summary(f2)$coefficients$cond[,"Estimate"],
    summary(f2)$coefficients$cond["age","Pr(>|z|)"]))                                    
  if('try-error' %in% class(res2)){
    res2 <- rep(NA,3)
  }
  names(res2) <- c("Beta0","Beta1","P-Beta")
  return(res2)
})

###############################################################
### define where the output is being written to 
outputTable.file <- file.path(my.workingDir, paste0("age-glmmTMB_",paramTest ,"_",cellTypeTest, "_resultsTable.txt"))
if(file.exists(outputTable.file)){
  write.table(t(res.t), file=outputTable.file, sep="\t", append=TRUE, row.names=TRUE, col.names=FALSE)
} else {
  write.table(t(res.t), file=outputTable.file, sep="\t", append=FALSE, row.names=TRUE, col.names=TRUE)
}