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
my.workingDir <- "/home/rd796/palmer_scratch/senage"

setwd(my.workingDir)

#### read in the arguments for the Rscript
inputArguments <- commandArgs(trailingOnly=TRUE)
cellTypeTest <- inputArguments[1]
genes.test<- inputArguments[2:length(inputArguments)]

####### load the seurat object
l=load('/home/rd796/project/ageproj/imm12_10.RData')

DefaultAssay(immune.combined)<-'RNA'

immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='Alv. Fibroblast']='Alv_Fibroblast'
immune.combined$predicted.id[immune.combined$predicted.id=='Alv. Macrophage']='Alv_Macrophage'
immune.combined$predicted.id[immune.combined$predicted.id=='Adventitial Fibroblast']='Adv_Fibroblast'

comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)
immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)

soup.subset <- subset(immune.combined, subset=predicted.id==cellTypeTest)

gc()

########### defining the parameters for yunqing's code
gene <- genes.test
subject <- soup.subset@meta.data$orig.ident
sf <- soup.subset@meta.data$nCount_RNA
senmayo <- soup.subset@meta.data$SenMayo1

### subset from the raw count matrix only the genes we are testing
count_mat <- soup.subset@assays$RNA@counts[gene, ]

### print to screen/sh.o what this job is gonna do
cat("Running GLMM\n", sep="")

############################################################### Ruben's code
res.t <- apply(count_mat, 1, function(gene){
  tmp.df <- data.frame(y=gene, norm_sf=sf, sm=senmayo, sub=as.numeric(factor(subject)))
  #f2 <- try(glmer(y ~ age + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = poisson))
  f2 <- try(glmmTMB(y ~ sm + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = nbinom2)) #zi = ~ age;
    res2 <-try(c(
    summary(f2)$coefficients$cond[,"Estimate"],
    summary(f2)$coefficients$cond["sm","Pr(>|z|)"]))                                    
  if('try-error' %in% class(res2)){
    res2 <- rep(NA,3)
  }
  names(res2) <- c("Beta0","Beta1","P-Beta")
  return(res2)
})

###############################################################
### define where the output is being written to 
outputTable.file <- file.path(my.workingDir, paste0("age-glmmTMB_", cellTypeTest, "_resultsTable.txt"))
if(file.exists(outputTable.file)){
  write.table(t(res.t), file=outputTable.file, sep="\t", append=TRUE, row.names=TRUE, col.names=FALSE)
} else {
  write.table(t(res.t), file=outputTable.file, sep="\t", append=FALSE, row.names=TRUE, col.names=TRUE)
}