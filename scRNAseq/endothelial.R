setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

#ENDOTHELIUM
Idents(immune.combined)=immune.combined@meta.data$predicted.id
endo <- subset(immune.combined, idents=comp[['Endothelial']])

DefaultAssay(endo)<-"integrated"
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 2000)
endo <- ScaleData(endo, verbose = FALSE)
endo <- RunPCA(endo, npcs = 30, verbose = FALSE)
endo <- RunUMAP(endo, reduction = "pca", dims = 1:20)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:20)
endo <- FindClusters(endo, resolution = 0.5)
DimPlot(endo,label=T,raster=F)+DimPlot(endo,group.by='predicted.id',raster=F,label=T)

#old plots
d0=DimPlot(endo,raster=F,label=T)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d1=DimPlot(endo,group.by='predicted.id',raster=F,label=T)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d2=DimPlot(endo,group.by='agebin',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d3=DimPlot(endo,group.by='orig.ident',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')+NoLegend()
d0+d1+d2+d3

#new plots
endo_df=fill_df[fill_df$comp=='Endothelial',]
p1<-DimPlot(endo,label=T,repel=T,group.by='predicted.id',cols=endo_df$color,order=rev(endo_df$predicted.id),shuffle=T)+NoLegend()+
        theme(strip.background = element_rect(fill="red"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p1[[1]][["data"]][["tempvar"]]='Cell Type'
p2<-DimPlot(endo,group.by='Manuscript_Identity',cols=c('#3E4E50','#FACFAD'))+
        theme(legend.position=c(.65,.75),
            strip.background = element_rect(fill="green"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title =element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p2[[1]][["data"]][["tempvar"]]='Dataset'
p3<-DimPlot(endo,group.by = 'agebin',cols=c('#91C4F2','#253C78'),shuffle=T)+
        theme(legend.position=c(.65,.75),
            strip.background = element_rect(fill="blue"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p3[[1]][["data"]][["tempvar"]]='Age'
plot_grid(p1+facet_wrap(~tempvar), p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar),ncol=1)

#SECOND ITERATION, GCAP ONLY
Idents(immune.combined)=immune.combined@meta.data$predicted.id
endo <- subset(immune.combined, idents='gCap')
DefaultAssay(endo)<-"integrated"
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 2000)
endo <- ScaleData(endo, verbose = FALSE)
endo <- RunPCA(endo, npcs = 30, verbose = FALSE)
endo <- RunUMAP(endo, reduction = "pca", dims = 1:20)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:20)
endo <- FindClusters(endo, resolution = 0.5)
DimPlot(endo,label=T,raster=F)+DimPlot(endo,group.by='agebin',raster=F,label=T)+DimPlot(endo,group.by='orig.ident',raster=F,label=F)
