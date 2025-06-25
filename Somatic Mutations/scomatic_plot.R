setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(ggpubr)

#LOAD DATA
load('agingseurat.RData')
immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
load('fill_df.RData')

##########
# GLOBAL #
##########
#Read glmmTMB output
datal=read.table(paste0('/home/rd796/palmer_scratch/global/age-glmmTMB_global_newmut2_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)
#datal=read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_GMUT_GENES/age-glmmTMB_global_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)
datal=datal[!(is.na(datal$Beta1)|is.na(datal$p_val)|!is.na(datal$blank)),]
datal=mutate(datal,p_val_adj=p.adjust(datal$p_val,method='fdr',n=nrow(datal)))
datal=datal[order(datal$p_val,decreasing=F),]
datal=datal[datal$p_val<0.05,]
noquote(datal$Gene[datal$Beta1>0][1:100])
noquote(datal$Gene[datal$Beta1<0&datal$p_val<0.05][1:100])

#Prepare for plotting
load('mutmatv2.RData')
ann <- gfg_data
ann$bin=factor(cut(ann$burden,c(0,median(ann$burden),100),c('Low Mutation','High Mutation')),levels=c('Low Mutation','High Mutation')) #OPTIONAL; BINARY MUT
small=subset(immune.combined,subset=orig.ident%in%rownames(ann))
small=AddMetaData(small,metadata=ann$mut[match(small$orig.ident,rownames(ann))],col.name='mut')

#CDKN2A Plot
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features='CDKN2A'))
datalist=data.frame(exp=t(avg),mut=ann$bin)
ggplot(datalist, aes(x=mut, y=exp))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+theme_classic()+theme(text = element_text(size=15))+
labs(x=NULL, y = "CDKN2A Expression")+scale_y_continuous(limits=c(0,0.35))+stat_compare_means(label='p.signif')

feats=c(datal$Gene[datal$Beta1>0][1:100],datal$Gene[datal$Beta1<0][1:100])
rowsplit=data.frame(Gene=rep(c('Increased','Decreased'),each=length(feats)/2))

small=subset(small,subset=orig.ident=='448',invert=T)
ann=ann[-which(rownames(ann)=='448'),]
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features=feats))

#Make Heatmap
#ann=ann[(rownames(ann) %in% substring(colnames(avg),5)),,F] #sanity
hanno=data.frame(Mutations=ann$burden, Age=ann$ages)
hcols=list(Age=circlize::colorRamp2(c(0, 100), c("white", "blue")))
colAnn <- HeatmapAnnotation(df = hanno, col=hcols,which='col',gap = unit(1, 'mm'))

rowfeats=c('SEL1L', 'SMAD5', 'ANAPC1', 'UBE4A', 'USP25', 'USP33', 'RAD50', 'MTREX', 'IFI16', 'SAMHD1', 'CBL', 'APC')
rowfeats=c(rowfeats,c('ETFB', 'COQ10A', 'MDH2', 'SCO2', 'IDH3B', 'COX4I1', 'GFER'))
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avg) %in% rowfeats), labels = rownames(avg)[rownames(avg)%in%rowfeats]))

avg2=t(apply(avg,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht=Heatmap(na.omit(avg2), name = "Expression",  column_split=data.frame(factor(hanno$'Mutations')),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F,cluster_row_slices = F, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = T,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=F,column_title = NULL, row_split=rowsplit,use_raster = F)
draw(ht, padding = unit(c(2, 2, 8, 4), "mm"),merge_legend=T)


################
# BY CELL TYPE #
################
#Read glmmTMB output
cell.types<-c('AT2','gCap')
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/cell/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_MUT_GENES/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})
names(datal)<-cell.types

#Prepare for plotting
load('cellmutmatv2.RData')
gfg_data=gfg_data[-which(rownames(gfg_data)=='424'),]
ht=NULL
for (i in cell.types){
ann <- gfg_data[,i]
rnames=rownames(gfg_data)
tmp=!is.na(ann)
#if(i=='gCap'){tmp[18]=F}
ann<-ann[tmp]
ann=factor(cut(ann,c(0,median(ann),100),c('Low Mutation','High Mutation')),levels=c('Low Mutation','High Mutation'))
rnames<-rnames[tmp]
ann=data.frame('Mutation Burden'=(ann)) #log10
rownames(ann)=rnames
#if(i=='gCap'){ann$Mutation.Burden[rownames(ann)=="460"]="High Mutation"}
ann=ann %>% mutate(Age=immune.combined$age[match(rownames(ann),immune.combined$orig.ident)],'Cell Type'=i)

small<-subset(immune.combined,subset=predicted.id==i&orig.ident %in% rnames)
small=AddMetaData(small,metadata=ann$'Mutation.Burden'[match(small$orig.ident,rnames)],col.name='mut')

VlnPlot(small,features='MT-CO1',split.by='orig.ident',group.by='mut')

feats=c(datal[[i]]$Gene[datal[[i]]$Beta1>0][1:100],datal[[i]]$Gene[datal[[i]]$Beta1<0][1:100])
rowsplit=data.frame(Gene=rep(c('Increased','Decreased'),each=length(feats)/2))
#feats=datal[[i]]$Gene[1:100]
#rowsplit=cut(sign(datal[[i]]$Beta1[1:100]),c(-10,0,10),c('Positive','Negative'))
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features=feats))

#Make Heatmap
ann=ann[(rownames(ann) %in% substring(colnames(avg),5)),,F] #sanity
names(ann)[1]<-'Mutations'
hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=circlize::colorRamp2(c(0, 100), c("white", "blue")),
'Mutations'=c('High Mutation'='#A20021','Low Mutation'='#F52F57'))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann, col=hcols,which='col',gap = unit(1, 'mm'))

noquote(head(datal[['gCap']]$Gene[datal[['gCap']]$Beta1>0],100))
rowfeats=c('USP24','CUL1','CUL3','CBL','SHPRH', 'EMSY','CREBBP')
rowfeats=c(rowfeats,c('VIPR1', 'FCN3', 'FENDRR', 'MT-ND2','MT-ND1','NDUFA1','NDUFA3','ATP5MG','ATP5F1E','ALDH2'))
rowfeats=c(rowfeats,c('INO80', 'PSIP1', 'METTL3', 'NOC3L'))
rowfeats=c(rowfeats,c('MACROD1', 'ABRAXAS2', 'BCCIP', 'RFC1', 'RAD17'))
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avg) %in% rowfeats), labels = rownames(avg)[rownames(avg)%in%rowfeats]))

avg2=t(apply(avg,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht0=Heatmap(na.omit(avg2), name = "Expression",  column_split=data.frame(factor(ann$'Mutations')),cluster_columns = T, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = TRUE,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=F,column_title = NULL, row_split=rowsplit,use_raster = F)
            ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 8, 6), "mm"),auto_adjust=F,merge_legend=T)
