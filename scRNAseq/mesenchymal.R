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
mes <- subset(immune.combined, idents=comp[['Mesenchymal']])

DefaultAssay(mes)<-"integrated"
mes <- FindVariableFeatures(mes, selection.method = "vst", nfeatures = 2000)
mes <- ScaleData(mes, verbose = FALSE)
mes <- RunPCA(mes, npcs = 30, verbose = FALSE)
mes <- RunUMAP(mes, reduction = "pca", dims = 1:20)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:20)
mes <- FindClusters(mes, resolution = 0.5)
DimPlot(mes,label=T,raster=F)+DimPlot(mes,group.by='predicted.id',raster=F,label=T)

#old plots
d1=DimPlot(mes,group.by='predicted.id',raster=F,label=T)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d2=DimPlot(mes,group.by='agebin',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')
d3=DimPlot(mes,group.by='orig.ident',raster=F,label=F)+ggtitle(NULL)+xlab('UMAP 1')+ylab('UMAP 2')+NoLegend()
d1+d2+d3

#new plots
mes_df=fill_df[fill_df$comp=='Mesenchymal',]
p1<-DimPlot(mes,raster=F,label=T,repel=T,group.by='predicted.id',cols=mes_df$color,order=rev(mes_df$predicted.id))+NoLegend()+
        theme(strip.background = element_rect(fill="red"),
              strip.text = element_text(size=15, colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p1[[1]][["data"]][["tempvar"]]='Cell Type'
p2<-DimPlot(mes,group.by='Manuscript_Identity',raster=F)+
        theme(legend.position=c(.65,.85),
            strip.background = element_rect(fill="green"),
              strip.text = element_text(size=15, colour="white"),
              plot.title =element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p2[[1]][["data"]][["tempvar"]]='Dataset'
p3<-FeaturePlot(mes,features = 'age',raster=F)+
        theme(legend.position=c(.85,.85),
            strip.background = element_rect(fill="blue"),
              strip.text = element_text(size=15, colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p3[[1]][["data"]][["tempvar"]]='Age'
plot_grid(p1+facet_wrap(~tempvar), p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar),ncol=3)
