setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

load('imm12_10.RData') #load integrated

#immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
#immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='Alv. Fibroblast']<-'Fibroblast'
immune.combined$predicted.id[immune.combined$predicted.id=='Adventitial Fibroblast']<-'Fibroblast'

comp<-list(Epithelial=c('AT1','AT2','AT2B','AT2S','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)

immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)
immune.combined<-subset(immune.combined,subset=predicted.id %in% c('AT2S','AT2B'),invert=F)

immune.combined<-subset(immune.combined,subset=Manuscript_Identity=='Baylor') #Or IPF Cell Atlas

output<-as.data.frame(immune.combined$predicted.id)
output2<-data.frame(Index=substring(rownames(output),0,16),Cell_type=output$`immune.combined$predicted.id`)
#output2<-data.frame(Index=substring(rownames(output),6),Cell_type=output$`immune.combined$predicted.id`) #IPF
write.table(output2, file='meta_scomatic.tsv', quote=FALSE, sep='\t', col.names = NA)
