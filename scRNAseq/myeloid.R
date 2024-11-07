setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

#Mesenchyme
Idents(immune.combined)=immune.combined@meta.data$predicted.id
mye <- subset(immune.combined, idents=comp[['Myeloid']])

DefaultAssay(mye)<-"integrated"
mye <- FindVariableFeatures(mye, selection.method = "vst", nfeatures = 2000)
mye <- ScaleData(mye, verbose = FALSE)
mye <- RunPCA(mye, npcs = 30, verbose = FALSE)
mye <- RunUMAP(mye, reduction = "pca", dims = 1:20)
mye <- FindNeighbors(mye, reduction = "pca", dims = 1:20)
mye <- FindClusters(mye, resolution = 0.5)
DimPlot(mye,label=T,raster=F)+DimPlot(mye,group.by='predicted.id',raster=F,label=T)
mye<-subset(mye,subset=integrated_snn_res.0.5==11,invert=T)
mye$predicted.id[Idents(mye) %in% c(3,9)]='Alv. Macrophage'

#old plots
d1=DimPlot(mye,group.by='predicted.id',raster=F,label=T)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d2=DimPlot(mye,group.by='agebin',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d3=DimPlot(mye,group.by='orig.ident',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')+NoLegend()
d1+d2+d3

#new plots
mye_df=fill_df[fill_df$comp=='Myeloid',]
p1<-DimPlot(mye,raster=F,label=T,repel=T,group.by='predicted.id',cols=rev(mye_df$color),order=(mye_df$predicted.id))+NoLegend()+
        theme(strip.background = element_rect(fill="red"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p1[[1]][["data"]][["tempvar"]]='Cell Type'
p2<-DimPlot(mye,group.by='orig.ident',raster=F)+
        theme(legend.position=c(.85,.1),
            strip.background = element_rect(fill="green"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title =element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())+NoLegend()
p2[[1]][["data"]][["tempvar"]]='Sample'
p3<-FeaturePlot(mye,features = 'age',raster=F)+
        theme(legend.position=c(-0,.12),
            strip.background = element_rect(fill="blue"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p3[[1]][["data"]][["tempvar"]]='Age'
plot_grid(p1+facet_wrap(~tempvar), p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar),ncol=3)
