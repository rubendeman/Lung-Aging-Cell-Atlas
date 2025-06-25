setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(stringr)

#read glmmTMB output
comp<-list(Epithelial=c('AT1','AT2B','AT2S','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-as.list(unlist(comp)); names(cell.types)<-cell.types
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/GLMM_nozi/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_AGE_GENES/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(datal,FUN = setNames,nm=c('Beta0','Beta1','p_val','p_val_adj'))
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})
names(datal)<-str_replace_all(names(datal), c("Alv_Macrophage" = "Alv. Macrophage", "Adv_Fibroblast" = "Adventitial Fibroblast", "Alv_Fibroblast"="Alv. Fibroblast"))
datalsig<-lapply(datal,function(x){x[x$p_val_adj<0.05,]})

#Load libra out
#load('/home/rd796/palmer_scratch/libra_out.RData')
load('/home/rd796/palmer_scratch/libra_out_hisample_at2sub.RData')
res.t<-res.t[order(res.t$p_val_adj),]
res.t<-res.t[res.t$p_val<1,]
datal2<-split(res.t,res.t$cell_type)
datal2sig<-lapply(datal2,function(x){x[x$p_val<0.05,]})

#Compare
cell.types=c('AT1','gCap','Aerocyte','Alv. Macrophage','SMC','Adventitial Fibroblast','Myofibroblast','Venous','Arterial','Lymphatic')
int=lapply(cell.types,function(x){intersect(datalsig[[x]]$Gene,datal2sig[[x]]$gene)})
names(int)<-cell.types
intsign1<-lapply(cell.types,function(x){sign(datalsig[[x]]$Beta1[match(int[[x]],datalsig[[x]]$Gene)])})
intsign2<-lapply(cell.types,function(x){sign(datal2sig[[x]]$avg_logFC[match(int[[x]],datal2sig[[x]]$gene)])})
names(intsign1)<-cell.types
names(intsign2)<-cell.types

unlist(lapply(cell.types,function(x){paste(x,nrow(datalsig[[x]]),length(int[[x]]))}))

#FET
for (x in cell.types){
f1=length(intersect(datalsig[[x]]$Gene,datal2sig[[x]]$gene)) #overlap
f2=length(intersect(datal[[x]]$Gene,datal2sig[[x]]$gene))-f1 #in a not b
f3=length(intersect(datalsig[[x]]$Gene,datal2[[x]]$gene))-f1 #in b not a
f4=length(intersect(datal[[x]]$Gene,datal2[[x]]$gene))-f1-f2-f3 #in b not a
matrix <- matrix(c(f1, f2, f3, f4), nrow=2)
print(fisher.test(matrix, alternative="greater"))
}

#UPSETS
library(UpSetR)
ct=comp[['Endothelial']]
u1<-upset(fromList(int[ct]),order.by='freq',empty.intersections="on")
u1

noquote(int[['SMC']][intsign1[['SMC']]!=1]) #this is positive genes
noquote(datalsig[['Myofibroblast']]$Gene[datalsig[['Myofibroblast']]$Beta1<0]) #this is positive genes
intersect(datalsig[['gCap']]$Gene[datalsig[['gCap']]$Beta1<0],datalsig[['Aerocyte']]$Gene[datalsig[['Aerocyte']]$Beta1<0])


#Load Seurat
load('agingseurat.RData')
#immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
#immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
load('fill_df.RData')



#BIG MAP BY SAMPLE (NEW)
cell.types='gCap'
mincells=3
ht=NULL
for (i in cell.types){
#feats=int[[i]]
#rowsplit=cut(intsign1[[i]],c(-10,0,10),c('Positive','Negative'))
feats=datalsig[[i]]$Gene
rowsplit=cut(sign(datalsig[[i]]$Beta1),c(-10,0,10),c('Positive','Negative'))
selgenes=list()
selgenes[['SMC']]=c('ITGA4','ITGA7','FGF1','COL4A6','LTBP1','LTBP2','ZNF44','ZNF83','ZNF138','ZNF148','ZNF423','ZNF444','CTBP2','SP1','GLI2','MECP2','MAML3')
selgenes[['Myofibroblast']]=c('TNFRSF1A','DUSP22','CARD19','C1S','PARP10','UBASH3B','DNASE2','FAM135A')
selgenes[['gCap']]=c('ABLIM1', 'MIB1', 'RNF19A', 'UBXN7', 'UBR3', 'UBE2H', 'UBE2K', 'UBE2E3', 'FBXL17', 'FBXO11', 'FBXO42', 'NPEPPS', 'VMP1', 'ULK1', 'RICTOR',
'SOS1', 'SOS2', 'MAP3K20', 'MAP4K4', 'MAPK8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND5', 'MT-CYB', 'VIPR1', 'IL7R', 'FCN3', 'HPGD')
selgenes[['Aerocyte']]=selgenes[['gCap']]
hl1 = which(feats %in% selgenes[[cell.types]])
hl2=1:10
ha = rowAnnotation(foo = anno_mark(at = c(hl1,hl2), labels = feats[c(hl1,hl2)]))
    noquote(feats)
#Downsample
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i))#&Manuscript_Identity=='Dataset 1')) #BAYLOR ONLY

#Hisamples only
cbmc.small$temp<-paste0(cbmc.small$predicted.id,"_",cbmc.small$orig.ident)
tbl<-as.data.frame(table(cbmc.small$temp))
keep<-tbl$Var1[tbl$Freq>mincells]
cbmc.small<-subset(cbmc.small,subset=temp %in% keep)

cbmc.small<-NormalizeData(cbmc.small)
Idents(cbmc.small)<-paste(cbmc.small$predicted.id,cbmc.small$agebin)
avgexp<-AverageExpression(cbmc.small,assays='RNA',layer='data',feats,group.by = c('predicted.id','orig.ident','agebin'))
avgexp<-as.matrix(avgexp$RNA)

#Scale Matrix
allind=as.data.frame(table(paste(cbmc.small$orig.ident,cbmc.small$Manuscript_Identity)))
allind=data.frame(do.call('rbind', strsplit(as.character(allind$Var1),' ',fixed=TRUE)))
ind1=na.omit(match(allind$X1[allind$X3==1],stringr::str_split_fixed(colnames(avgexp),"_",2)[,1]))
ind2=na.omit(match(allind$X1[allind$X3==2],stringr::str_split_fixed(colnames(avgexp),"_",2)[,1]))
#avgexp[,ind1]=t(apply(avgexp[,ind1],1,function(x){(x-min(x))/(max(x)-min(x))}))
#avgexp[,ind2]=t(apply(avgexp[,ind2],1,function(x){(x-min(x))/(max(x)-min(x))}))
avgexp=cbind(t(apply(avgexp[,ind1],1,function(x){(x-min(x))/(max(x)-min(x))})),t(apply(avgexp[,ind2],1,function(x){(x-min(x))/(max(x)-min(x))})))
avgexp[is.na(avgexp)]<-0

#Create Annotations
    ann <- data.frame(i,sub('_.*','',colnames(avgexp)),sub('.*_','',colnames(avgexp)))
    colnames(ann) <- c('Cell Type','Sample','Age')
    ann<-ann %>% mutate(Dataset=cbmc.small$Manuscript_Identity[match(ann$Sample,cbmc.small$orig.ident)])
    ann<-ann[,c(1,3,4,2)]
#ann <- data.frame(sub('_.*','',colnames(avgexp)),cut(cbmc.small$age[match(gsub('.*_(.+)_.*','\\1',colnames(avgexp)),cbmc.small$orig.ident)],c(0,40,60,100),c('young','mid','old')))

hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))

col_fun = circlize::colorRamp2(c(0.2, 1), c("black", "#FFFF00"))
#col_fun = circlize::colorRamp2(c(0,0.5, 1), c("purple", "black", "#FFFF00"))
ht0=Heatmap(na.omit(avgexp), name = "Expression",  column_split=data.frame(factor(ann$'Age',levels=c('Young','Aged')),factor(ann$Dataset)),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
    cluster_row_slices = F, cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = T,
    show_row_dend = F, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
    show_column_names = F, show_row_names=F, row_split=rowsplit, use_raster = F)
    ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 2, 5), "mm"),merge_legend=T,ht_gap=unit(0.5,'cm'))
