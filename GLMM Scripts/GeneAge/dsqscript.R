.libPaths(c("/home/rd796/project/R/4.2", .libPaths()))
library(Seurat)
library(tidyverse)

####### set working directory
my.workingDir <- "/home/rd796/palmer_scratch"

setwd(my.workingDir)

load('/home/rd796/project/ageproj/imm12_10.RData')
DefaultAssay(immune.combined)<-'RNA'
genes.test <- rownames(immune.combined)[rowSums(immune.combined@assays$RNA@counts>0)>100] #keep genes expressed in >100 cells

#### we need to split the genes up across jobs

gene.groups <- split(genes.test, ceiling(seq_along(genes.test)/100))
# all(unlist(gene.groups) %in% genes.test)
# all(genes.test %in% unlist(gene.groups))

#### define where all the scripts are going to be written to
script.output.dir <- file.path(my.workingDir,"glmm.scripts")
### make this folder if it doesn't already exist
if(dir.exists(script.output.dir)==FALSE){
  dir.create(script.output.dir)
}

jobsub.filepath <- file.path(script.output.dir, "joblist.txt")

rscript.function.filepath <- "/gpfs/gibbs/pi/kaminski/public/Backup/Ruben/GLMMAging/runmodel.R"

### what cell type to test
comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)
#cellType=cell.types
cellType=c('AT2','gCap')

### create a new jobsub file before starting
cat("", sep="", file=jobsub.filepath, append=FALSE)

cmd.out <- NULL
for(i in 1:length(gene.groups)){
  for(j in 1:length(cellType)){
  temp.cellType <- cellType[j]
  temp.genes <- gene.groups[[i]] 
  cmd.out <- paste("module load R/4.2.0-foss-2020b; ", sep="")
  cmd.out <- paste(cmd.out, "Rscript ", rscript.function.filepath, " ", paste(c(temp.cellType, temp.genes), collapse=" "), sep="")
  cat(cmd.out, "\n", sep="", file=jobsub.filepath, append=TRUE)
  }
}
