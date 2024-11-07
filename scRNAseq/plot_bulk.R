setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(ComplexHeatmap)
library(stringr)
library(WGCNA)

#LOAD DATA
load('loggedmes.RData') #load module eigenegenes
load('prop12_15.RData') #load music props
comp<-list(epi=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),endo=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),mes=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),imm=c('Monocyte','Macrophage','Alv. Macrophage'),lymph=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)

proportion.df2<-tempprop$Est.prop.weighted
proportion.df2=proportion.df2[,cell.types]

#Bulk
TraitCor = cor(MEs, proportion.df2, use = "p")
TraitPvalue = corPvalueStudent(TraitCor, 572)
menames=read.csv('menames.csv')
mtc2=rownames(TraitCor)
menames=as.matrix(menames)
#menames=c('ME0',menames)
mtc=match(menames,mtc2)
SortCor=TraitCor[mtc,]
rownames(SortCor)<-substring(rownames(SortCor),3)

SortCor[is.na(SortCor)]<-0
quantile(SortCor, c(0.1, 0.9))

#scRNAseq
load('geneInfo0.RData') #load gene modules
load('agecorgenes.RData') #load genes
datalist2 = list()
DefaultAssay(immune.combined)<-"RNA"
for (i in agecorgenes){ #from wgcna_output.R
#1:ncol(MEs)-1
#temp=geneInfo0$geneSymbol[geneInfo0$moduleColor==i]
#immune.combined <- ScaleData(object = immune.combined, features = as.matrix(temp))
#geneList <- list(as.matrix(temp))
#immune.combined <- AddModuleScore(immune.combined, features = geneList, name=as.character(i))
m=data.frame(immune.combined@meta.data$predicted.id,immune.combined@meta.data$age,dplyr::select(immune.combined@meta.data,paste0("X",i,"1")))
o=m %>% group_by(immune.combined.meta.data.predicted.id) %>% summarise(V1=mean(get(paste0("X",i,"1"))))
datalist2[[i]]=o$V1
}
DefaultAssay(immune.combined)<-"integrated"
datalist2=do.call("cbind",datalist2)
rownames(datalist2)=o$immune.combined.meta.data.predicted.id
colnames(datalist2)=paste(agecorgenes)

listage2=((scale(datalist2)[cell.types,]))#[,agecorgenes]
listsm2=((scale(datalist2[cell.types,c('2','15','19','20','28','31')])))

#BOTH
i1=SortCor[agecorgenes,] #from wgcna_output.R
h=Heatmap(t(listage2), name = "Expression", 
          cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = T,
          show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00")))ord<-row_order(h)
h2=Heatmap((i1), name = "Correlation", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = T, 
           show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#FF00FF", "black", "#FFFF00")))
draw(h %v% h2, padding = unit(c(10, 2, 2, 10), "mm"))

h=Heatmap(t(listage2), name = "Expression", 
          cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = T,
          show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00")))ord<-row_order(h)
h2=Heatmap((i1), name = "Correlation", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = T, 
           show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#FF00FF", "black", "#FFFF00")))
draw(h + h2, padding = unit(c(10, 2, 2, 10), "mm"))

#sen
i2=SortCor[c(2,15,19,20,28,31),]
h=Heatmap(t(listsm2), name = "Expression", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = FALSE,
         show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00")))
h2=Heatmap((i2), name = "Correlation", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = FALSE,
         show_row_dend = FALSE, row_names_gp = gpar(fontsize = 8), col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#FF00FF", "black", "#FFFF00")))
draw(h %v% h2, padding = unit(c(8, 2, 2, 8), "mm"))
#or invert
i2=SortCor[c(2,15,19,20,28,31),]
h=Heatmap((listsm2), name = "Expression", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = FALSE, row_names_side='left',
         show_row_dend = FALSE, row_names_gp = gpar(fontsize = 12), col = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00")))
h2=Heatmap(t(i2), name = "Correlation", cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = FALSE, row_names_side='left',
         show_row_dend = FALSE, row_names_gp = gpar(fontsize = 12), col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#FF00FF", "black", "#FFFF00")))
draw(h + h2, padding = unit(c(2, 8, 2, 8), "mm"))
