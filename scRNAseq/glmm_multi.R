setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage','DC'),Lymphoid=c('B','T','Mast','NK'))
cell.types<-unlist(comp)

#read glmmTMB output
cell.types<-c('AT1','AT2','gCap')
comparisons<-c('Late-Aged','Early-Aged','Male','Female','1','2','Yes','No')
datal=list()
for(k in comparisons){
for(j in cell.types){
    datal[[paste(j,k)]]<-read.table(paste0('/home/rd796/palmer_scratch/multi/age-glmmTMB_',k,'_',j,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)
}
}
#datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_AGE_GENES/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(datal,FUN = setNames,nm=c('Beta0','Beta1','p_val','p_val_adj'))
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})

datasig<-lapply(datal,function(x){x[x$p_val<0.05,]})
datafdr<-lapply(datal,function(x){x[x$p_val_adj<0.05,]})

datasig10<-lapply(datal,function(x){x$Gene[x$p_val<0.05&x$Beta1>0][1:10]})

#Compare
int=list()
intsign1=list()
intsign2=list()
comparisons<-list(agetri=c('Late-Aged','Early-Aged'),sex=c('Male','Female'),Manuscript_Identity=c('1','2'),smoke=c('Yes','No'))
for(j in cell.types){
for(cond in names(comparisons)){
    k=comparisons[[cond]]
int[[paste(j,cond)]]<-intersect(datasig[[paste(j,k[1])]]$Gene,datasig[[paste(j,k[2])]]$Gene)
intsign1[[paste(j,cond)]]<-sign(as.numeric(datasig[[paste(j,k[1])]]$Beta1[match(int[[paste(j,cond)]],datasig[[paste(j,k[1])]]$Gene)]))
intsign2[[paste(j,cond)]]<-sign(as.numeric(datasig[[paste(j,k[2])]]$Beta1[match(int[[paste(j,cond)]],datasig[[paste(j,k[2])]]$Gene)]))
}
}

#Unique
uni<-lapply(1:length(datasig),function(x){datasig[[x]][!datasig[[x]]$Gene %in% int[[ceiling(x/2)]],]})
names(uni)<-names(datasig)
unifdr<-lapply(uni,function(x){x[x$p_val_adj<0.05&x$Beta1>0,]})
uni10<-lapply(uni,function(x){x$Gene[x$p_val<0.05&x$Beta1>0][1:30]})

noquote(unifdr[['gCap Female']]$Gene) #this is positive genes


#LOAD DATA
load('agingseurat.RData')
immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
immune.combined$agetri<-as.character(immune.combined$agetri)
immune.combined$agetri[immune.combined$agetri=='Middle-Aged']<-'Early-Aged'
immune.combined$agetri[immune.combined$agetri=='Aged']<-'Late-Aged'
immune.combined$agetri<-factor(immune.combined$agetri,levels=c('Young','Early-Aged','Late-Aged'))
load('fill_df.RData')



#BIG MAP BY SAMPLE (NEW)
i=1
cell.types=c('AT1','AT2','gCap')
freqs=c(3,10,10); names(freqs)=cell.types
cell.types=cell.types[i]
freqs=freqs[i]
var='agetri'
#remove duplicate genes
conds0=expand.grid(cell.types,comparisons[[var]])
conds<-paste(conds0$Var1,conds0$Var2)
picklist<-lapply(uni10[conds],data.frame)
picklist=bind_rows(picklist,.id='source')
colnames(picklist)<-c('source','gene')
picklist=picklist[match(unique(picklist$gene),picklist$gene),]
nrsplit0=picklist %>% group_by (source) %>% summarise(length(source))
colnames(nrsplit0)<-c('source','length')
nrsplit1=nrsplit0$length[match(conds,nrsplit0$source)]

ht=NULL
for (i in cell.types){
feats=picklist$gene
rowsplit=rep(conds,nrsplit1)
selgenes=list()
selgenes[['agetri']]=c('DNMT3A', 'BRCA2', 'TP53', 'IRF3', 'FBXL17', 'FBXO42', 'ABLIM1', 'MIB1', 'ANAPC1')
selgenes[['smoke']]=c('RAD21', 'INO80', 'CTDSP2', 'ZMYM3', 'PRKCD', 'CYBC1', 'PRKCD', 'IL17RE')
selgenes[['gCap']]=c('ABLIM1', 'MIB1', 'RNF19A', 'UBXN7', 'UBR3', 'UBE2H', 'UBE2K', 'UBE2E3', 'FBXL17', 'FBXO11', 'FBXO42', 'NPEPPS', 'VMP1', 'ULK1', 'RICTOR',
'SOS1', 'SOS2', 'MAP3K20', 'MAP4K4', 'MAPK8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND5', 'MT-CYB', 'VIPR1', 'IL7R', 'FCN3', 'HPGD')
load('senmayolist.RData')
selgenes[['senmayo']]=unlist(senmayolist)
hl1 = which(feats %in% unlist(selgenes))
ha = rowAnnotation(foo = anno_mark(at = c(hl1), labels = feats[c(hl1)]))
    noquote(feats)
#Downsample
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i))#&Manuscript_Identity=='Dataset 1')) #BAYLOR ONLY

#Hisamples only
cbmc.small$temp<-paste0(cbmc.small$predicted.id,"_",cbmc.small$orig.ident)
tbl<-as.data.frame(table(cbmc.small$temp))
keep<-tbl$Var1[tbl$Freq>freqs[i]]
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
    temp<-cbmc.small@meta.data[,var]
    ann<-ann %>% mutate(Var0=temp[match(ann$Sample,cbmc.small$orig.ident)])
    ann<-ann[,c(1,3,4,2)]
#ann <- data.frame(sub('_.*','',colnames(avgexp)),cut(cbmc.small$age[match(gsub('.*_(.+)_.*','\\1',colnames(avgexp)),cbmc.small$orig.ident)],c(0,40,60,100),c('young','mid','old')))

vec=c('purple','pink')
names(vec)<-comparisons[[var]]
if(var=='smoke'){vec=c(vec,Unknown="#ADCAD6")}
if(var=='agetri'){vec=c(vec,Young="#ADCAD6")}
hcols=list('Cell Type'=fill_df$color[match(cell.types,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'),Var0=vec)
names(hcols$'Cell Type')=cell.types
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))
ann2<-conds0[rep(seq_len(nrow(conds0)),nrsplit1),]
colnames(ann2)<-c('Cell Type','Var0')
rowAnn <- HeatmapAnnotation(df = ann2,col=hcols,which = 'row',gap = unit(1, 'mm'),show_legend=c(T,T))
if(i!=cell.types[1]){
    rowAnn=NULL
}
if(length(cell.types)>1){
if(i!=cell.types[3]){
    ha=NULL
}
}

col_fun = circlize::colorRamp2(c(0.2, 1), c("black", "#FFFF00"))
#col_fun = circlize::colorRamp2(c(0,0.5, 1), c("purple", "black", "#FFFF00"))
ht0=Heatmap(na.omit(avgexp), name = "Expression",  column_split=data.frame(ann$Var0,factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,left_annotation=rowAnn,right_annotation=ha,
    cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
    show_row_dend = F, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
    show_column_names = F, show_row_names=F, row_split=rowsplit, use_raster = F)
    ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 2, 5), "mm"),merge_legend=T,ht_gap=unit(0.1,'cm'))



#BIG MAP BY SAMPLE (NEWER)
cell.types=c('AT1','AT2','gCap')
freqs=c(3,10,10); names(freqs)=cell.types
var='agetri'
ht=NULL
for (i in cell.types){
    #remove duplicate genes
conds0=expand.grid(i,comparisons[[var]])
conds<-paste(conds0$Var1,conds0$Var2)
picklist<-lapply(uni10[conds],data.frame)
picklist=bind_rows(picklist,.id='source')
colnames(picklist)<-c('source','gene')
picklist=picklist[match(unique(picklist$gene),picklist$gene),]
nrsplit0=picklist %>% group_by (source) %>% summarise(length(source))
colnames(nrsplit0)<-c('source','length')
nrsplit1=nrsplit0$length[match(conds,nrsplit0$source)]


feats=picklist$gene
rowsplit=rep(conds,nrsplit1)
selgenes=list()
selgenes[['agetri']]=c('DNMT3A', 'BRCA2', 'TP53', 'IRF3', 'FBXL17', 'FBXO42', 'ABLIM1', 'MIB1', 'ANAPC1')
selgenes[['smoke']]=c('RAD21', 'INO80', 'CTDSP2', 'ZMYM3', 'PRKCD', 'CYBC1', 'PRKCD', 'IL17RE')
selgenes[['gCap']]=c('ABLIM1', 'MIB1', 'RNF19A', 'UBXN7', 'UBR3', 'UBE2H', 'UBE2K', 'UBE2E3', 'FBXL17', 'FBXO11', 'FBXO42', 'NPEPPS', 'VMP1', 'ULK1', 'RICTOR',
'SOS1', 'SOS2', 'MAP3K20', 'MAP4K4', 'MAPK8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND5', 'MT-CYB', 'VIPR1', 'IL7R', 'FCN3', 'HPGD')
load('senmayolist.RData')
selgenes[['senmayo']]=unlist(senmayolist)
hl1 = which(feats %in% unlist(selgenes))
hl2=which(feats %in% c('LRP10','DUSP7','MAP2K2','NEK1','TBL2','BTC','HSPA6','ANXA4','STAT5B'))
ha = rowAnnotation(foo = anno_mark(at = c(hl1,hl2), labels = feats[c(hl1,hl2)]))
    noquote(feats)
#Downsample
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i))#&Manuscript_Identity=='Dataset 1')) #BAYLOR ONLY

#Hisamples only
cbmc.small$temp<-paste0(cbmc.small$predicted.id,"_",cbmc.small$orig.ident)
tbl<-as.data.frame(table(cbmc.small$temp))
keep<-tbl$Var1[tbl$Freq>freqs[i]]
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
    temp<-cbmc.small@meta.data[,var]
    ann<-ann %>% mutate(Var0=temp[match(ann$Sample,cbmc.small$orig.ident)])
    ann<-ann[,c(1,3,4,2)]
#ann <- data.frame(sub('_.*','',colnames(avgexp)),cut(cbmc.small$age[match(gsub('.*_(.+)_.*','\\1',colnames(avgexp)),cbmc.small$orig.ident)],c(0,40,60,100),c('young','mid','old')))

vec=c('purple','pink')
names(vec)<-comparisons[[var]]
if(var=='smoke'){vec=c(vec,Unknown="#ADCAD6")}
if(var=='agetri'){vec=c(vec,Young="#ADCAD6")}
hcols=list('Cell Type'=fill_df$color[match(cell.types,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'),Var0=vec)
names(hcols$'Cell Type')=cell.types
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))
ann2<-conds0[rep(seq_len(nrow(conds0)),nrsplit1),]
colnames(ann2)<-c('Cell Type','Var0')
rowAnn <- HeatmapAnnotation(df = ann2,col=hcols,which = 'row',gap = unit(1, 'mm'),show_legend=c(T,T))

col_fun = circlize::colorRamp2(c(0.2, 1), c("black", "#FFFF00"))
#col_fun = circlize::colorRamp2(c(0,0.5, 1), c("purple", "black", "#FFFF00"))
ht0=Heatmap(na.omit(avgexp), name = "Expression",  column_split=data.frame(ann$Var0,factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,left_annotation=rowAnn,right_annotation=ha,
    cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
    show_row_dend = F, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
    show_column_names = F, show_row_names=F, row_split=rowsplit, row_title=NULL, use_raster = F)
    ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 2, 5), "mm"),merge_legend=T,ht_gap=unit(0.1,'cm'))
