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

#ADD METADATA
#age
subnames=soup.subset@meta.data$subject.ident
ages=read.csv('sc_demo.csv') #load demographic file
mtch=match(subnames,ages[,1])
ages2=ages[,3][mtch]
soup.subset <- AddMetaData(object = soup.subset, metadata = ages2, col.name = 'age')

#agebin
ages3  <- vector(mode='character',length=length(ages2))
for (x in 1:length(ages2)) {
    val=ages2[x]
    if (val<61){ages3[x]="young";} else {ages3[x]="old"
    }
}
soup.subset <- AddMetaData(object = soup.subset, metadata = ages3, col.name = 'agebin')

#PREPROCESSING
soup.subset <- NormalizeData(soup.subset)
soup.subset <- FindVariableFeatures(soup.subset, selection.method = "vst", nfeatures = 2000)
soup.subset <- ScaleData(soup.subset, verbose = FALSE)
soup.subset <- RunPCA(soup.subset, npcs = 30, verbose = FALSE)
soup.subset <- RunUMAP(soup.subset, reduction = "pca", dims = 1:20)
soup.subset <- FindNeighbors(soup.subset, reduction = "pca", dims = 1:20)
soup.subset <- FindClusters(soup.subset, resolution = 0.5)
Idents(soup.subset)=soup.subset@meta.data$Manuscript_Identity

library('stringr');
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Capillary_A","Aerocyte");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Capillary_B","gCap");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"cDC1","Dendritic");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"cMonocyte","Monocyte");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"Macrophage_Alveolar","Alv. Macrophage");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Arterial","Arterial");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Venous","Venous");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"VE_Peribronchial","Peribronchial");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"B_Plasma","Plasma B");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"T_Cytotoxic","Cytotoxic T");
soup.subset@meta.data$Manuscript_Identity=str_replace(soup.subset@meta.data$Manuscript_Identity,"T_Regulatory","Regulatory T");
Idents(soup.subset) <-  soup.subset@meta.data$Manuscript_Identity

#SENESCENCE SCORING
#sen
rr=read.table('listsenmayo.txt'); #mapping high senescence score cells
rrr=c(rr[1:18,],rr[20:125,]);
soup.subset <- ScaleData(object = soup.subset, features = rrr);
geneList <- list(rrr);
soup.subset <- AddModuleScore(soup.subset, features = geneList, name="SenMayo");

#senbin
range(soup.subset@meta.data$SenMayo1)
sum(soup.subset@meta.data$SenMayo1>0.24)
ages3  <- vector(mode='character',length=length(ages2))
for (x in 1:length(ages2)) {
    val=soup.subset@meta.data$SenMayo1[x]; if (val<0.22){ages3[x]="Low Senescence Score";} else {ages3[x]="High Senescence Score";}
}
table(ages3)
soup.subset <- AddMetaData(object = soup.subset, metadata = ages3, col.name = 'senbin')

#REMOVE CELL TYPES
Idents(soup.subset)=soup.subset@meta.data$Manuscript_Identity
all.keep <- (c('ATI','ATII','Basal','Ciliated','Fibroblast','Myofibroblast','Monocyte','Macrophage','Alv. Macrophage','B','T','Arterial','Venous','Aerocyte','gCap'))
soup.subset <- subset(soup.subset, idents=all.keep)
