setwd('/home/rd796/project/ageproj')
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(Matrix)

#PREPROCESSING
sids<-c('273','275','302','330','338','356','367','368','371','376','383','396','398','399','408','415','416','424','429','435','445','446','448','454','460','475','478','484','485','491','499','501')
dir1 = "BCM_controls/BCM_"
dir2 = "/filtered_feature_bc_matrix"
counter=0;
matlist=list()
for (i in sids){
    counter=counter+1;
l1=Read10X(paste0(dir1,i,dir2))
matlist[[counter]]<-CreateSeuratObject(counts = l1,project=i)
}
#lungaging <- merge(matlist[[1]], y = unlist(matlist)[2:32], add.cell.ids = sids)

#RECIPROCAL INTEGRATION
matlist <- lapply(X = matlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = matlist)
matlist <- lapply(X = matlist, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = matlist, reference = c(8, 21), reduction = "rpca", dims = 1:30)
bm40k.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

VlnPlot(bm40k.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(bm40k.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
bm40k.integrated <- subset(bm40k.integrated, subset = nFeature_RNA > 750 & nFeature_RNA < 7500)

bm40k.integrated <- ScaleData(bm40k.integrated, verbose = FALSE)
bm40k.integrated <- RunPCA(bm40k.integrated, npcs=30, verbose = FALSE) #SPECIFY #pcs
bm40k.integrated <- RunUMAP(bm40k.integrated, dims = 1:30)

#save(bm40k.integrated,file='baylor_int.RData')

l=load('baylor_int.RData')
lungaging <- bm40k.integrated

#CELL TYPE LABEL TRANSFER
load('ppsoup.RData') #load IPF cell atlas

cell.anchors <- FindTransferAnchors(reference = soup.subset, query = lungaging, reduction="rpca", dims = 1:30, k.filter=NA)
predictions <- TransferData(anchorset = cell.anchors, refdata = soup.subset$Manuscript_Identity, dims = 1:30)
lungaging <- AddMetaData(lungaging, metadata = predictions)

#INTEGRATION
ifnb.list=list(lungaging,soup.subset)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features,reduction="rpca",k.filter=NA)
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

#QUALITY CONTROL VERIFY
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, assay='RNA', pattern = "^MT-")
immune.combined <- subset(immune.combined, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 10)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by='Manuscript_Identity',pt.size=0)

#RUN UMAP
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
#immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
#immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#UPDATE METADATA
immune.combined@meta.data$predicted.id<-coalesce(immune.combined@meta.data$predicted.id,immune.combined@meta.data$Manuscript_Identity)

immune.combined@meta.data$orig.ident[immune.combined@meta.data$orig.ident=='CellSoup']<-NA
immune.combined@meta.data$orig.ident<-coalesce(immune.combined@meta.data$orig.ident,immune.combined@meta.data$subject.ident)

ind<-is.na(immune.combined@meta.data$subject.ident)
immune.combined@meta.data$Manuscript_Identity[ind]<-"Baylor"
immune.combined@meta.data$Manuscript_Identity[!ind]<-"IPF Cell Atlas"
