setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
load('soup.subset.Robj') #load the IPF cell atlas

#VIEW DATASET
infolist=soup.subset@meta.data %>%
     group_by(subject.ident, Manuscript_Identity) %>%
     summarize("counts"=n()) %>%
     spread(value=counts, key=Manuscript_Identity) %>%
     t() %>% as.data.frame()

#PREPROCESSING
soup.subset <- NormalizeData(soup.subset)
soup.subset <- FindVariableFeatures(soup.subset, selection.method = "vst", nfeatures = 2000)
soup.subset <- ScaleData(soup.subset, verbose = FALSE)
soup.subset <- RunPCA(soup.subset, npcs = 30, verbose = FALSE)
soup.subset <- RunUMAP(soup.subset, reduction = "pca", dims = 1:20)
soup.subset <- FindNeighbors(soup.subset, reduction = "pca", dims = 1:20)
soup.subset <- FindClusters(soup.subset, resolution = 0.5)
Idents(soup.subset)=soup.subset@meta.data$Manuscript_Identity

library('stringr')
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Capillary_A","Aerocyte")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Capillary_B","gCap")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"cDC1","Dendritic")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"cMonocyte","Monocyte")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"Macrophage_Alveolar","Alv. Macrophage")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Arterial","Arterial")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Venous","Venous")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Peribronchial","Peribronchial")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"B_Plasma","Plasma B")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"T_Cytotoxic","Cytotoxic T")
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"T_Regulatory","Regulatory T")
Idents(soup.subset) <-  soup.subset@meta.data$Manuscript_Identity

#REMOVE CELL TYPES
Idents(soup.subset)=soup.subset@meta.data$Manuscript_Identity
all.keep <- (c('ATI','ATII','Basal','Ciliated','Fibroblast','Myofibroblast','Monocyte','Macrophage','Alv. Macrophage','B','T','Arterial','Venous','Aerocyte','gCap'))
soup.subset <- subset(soup.subset, idents=all.keep)
