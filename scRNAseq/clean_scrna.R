setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

l=load('imm10_29.RData')

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.2)
DimPlot(immune.combined,label=T,raster=F)+DimPlot(immune.combined,group.by='predicted.id',raster=F,label=T)

#EPITHELIUM
epi <- subset(immune.combined, subset=integrated_snn_res.0.2 %in% c(4,7,11,12))
myeloid <- subset(immune.combined, subset=integrated_snn_res.0.2 %in% c(0,1,3,6,8))
lymphoid <- subset(immune.combined, subset=integrated_snn_res.0.2 %in% c(2,13)) 
endo <- subset(immune.combined, subset=integrated_snn_res.0.2 %in% c(5)) 
mes <- subset(immune.combined, subset=integrated_snn_res.0.2 %in% c(9,14)) 

immune.combined@meta.data$predicted.id[immune.combined@active.ident==15]<-'Lymphatic'
immune.combined@meta.data$predicted.id[immune.combined@active.ident==10]<-'Multiplet'

DefaultAssay(epi)<-"integrated"
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:20)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:20)
epi <- FindClusters(epi, resolution = 0.5)
DimPlot(epi,label=T,raster=F)+DimPlot(epi,group.by='predicted.id',raster=F,label=T)
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(1,4,7,8,18)],rownames(immune.combined@meta.data))]<-'Ciliated'
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(15)],rownames(immune.combined@meta.data))]<-'Basal'
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(11)],rownames(immune.combined@meta.data))]<-'AT2S'
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(9)],rownames(immune.combined@meta.data))]<-'AT2B'
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(5)],rownames(immune.combined@meta.data))]<-'AT1'
immune.combined@meta.data$predicted.id[match(rownames(epi@meta.data)[epi@active.ident %in% c(3,12,13,14,16,17)],rownames(immune.combined@meta.data))]<-'Multiplet'

DefaultAssay(myeloid)<-"integrated"
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)
myeloid <- ScaleData(myeloid, verbose = FALSE)
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.1)
DimPlot(myeloid,label=T,raster=F)+DimPlot(myeloid,group.by='predicted.id',raster=F,label=T)
immune.combined@meta.data$predicted.id[match(rownames(myeloid@meta.data)[myeloid@active.ident %in% c(1,4)],rownames(immune.combined@meta.data))]<-'Alv. Macrophage'
immune.combined@meta.data$predicted.id[match(rownames(myeloid@meta.data)[myeloid@active.ident %in% c(0,5)],rownames(immune.combined@meta.data))]<-'Macrophage'
immune.combined@meta.data$predicted.id[match(rownames(myeloid@meta.data)[myeloid@active.ident %in% c(2)],rownames(immune.combined@meta.data))]<-'Monocyte'

DefaultAssay(lymphoid)<-"integrated"
lymphoid <- FindVariableFeatures(lymphoid, selection.method = "vst", nfeatures = 2000)
lymphoid <- ScaleData(lymphoid, verbose = FALSE)
lymphoid <- RunPCA(lymphoid, npcs = 30, verbose = FALSE)
lymphoid <- RunUMAP(lymphoid, reduction = "pca", dims = 1:20)
lymphoid <- FindNeighbors(lymphoid, reduction = "pca", dims = 1:20)
lymphoid <- FindClusters(lymphoid, resolution = 0.3)
DimPlot(lymphoid,label=T,raster=F)+DimPlot(lymphoid,group.by='predicted.id',raster=F,label=T)
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(7)],rownames(immune.combined@meta.data))]<-'B'
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(5)],rownames(immune.combined@meta.data))]<-'DC'
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(8)],rownames(immune.combined@meta.data))]<-'Mast'
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(6)],rownames(immune.combined@meta.data))]<-'Multiplet'
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(3)],rownames(immune.combined@meta.data))]<-'NK'
immune.combined@meta.data$predicted.id[match(rownames(lymphoid@meta.data)[lymphoid@active.ident %in% c(0,1,2,4,9)],rownames(immune.combined@meta.data))]<-'T'

DefaultAssay(endo)<-"integrated"
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 2000)
endo <- ScaleData(endo, verbose = FALSE)
endo <- RunPCA(endo, npcs = 30, verbose = FALSE)
endo <- RunUMAP(endo, reduction = "pca", dims = 1:20)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:20)
endo <- FindClusters(endo, resolution = 0.5)
DimPlot(endo,label=T,raster=F)+DimPlot(endo,group.by='predicted.id',raster=F,label=T)
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(7,9)],rownames(immune.combined@meta.data))]<-'Multiplet'
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(0,2,4,5)],rownames(immune.combined@meta.data))]<-'gCap'
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(1)],rownames(immune.combined@meta.data))]<-'Aerocyte'
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(3)],rownames(immune.combined@meta.data))]<-'Arterial'
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(6)],rownames(immune.combined@meta.data))]<-'Venous'
immune.combined@meta.data$predicted.id[match(rownames(endo@meta.data)[endo@active.ident %in% c(8)],rownames(immune.combined@meta.data))]<-'Peribronchial'

DefaultAssay(mes)<-"integrated"
mes <- FindVariableFeatures(mes, selection.method = "vst", nfeatures = 2000)
mes <- ScaleData(mes, verbose = FALSE)
mes <- RunPCA(mes, npcs = 30, verbose = FALSE)
mes <- RunUMAP(mes, reduction = "pca", dims = 1:20)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:20)
mes <- FindClusters(mes, resolution = 0.1)
DimPlot(mes,label=T,raster=F)+DimPlot(mes,group.by='predicted.id',raster=F,label=T)
immune.combined@meta.data$predicted.id[match(rownames(mes@meta.data)[mes@active.ident %in% c(4,5,6)],rownames(immune.combined@meta.data))]<-'Multiplet'
