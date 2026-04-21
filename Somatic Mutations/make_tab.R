library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

#LOAD INTEGRATED OBJECT

immune.combined<-subset(immune.combined,subset=Manuscript_Identity=='Baylor') #Or IPF Cell Atlas

output<-as.data.frame(immune.combined$predicted.id)
output2<-data.frame(Index=substring(rownames(output),0,16),Cell_type=output$`immune.combined$predicted.id`)
#output2<-data.frame(Index=substring(rownames(output),6),Cell_type=output$`immune.combined$predicted.id`) #IPF
write.table(output2, file='meta_scomatic.tsv', quote=FALSE, sep='\t', col.names = NA)
