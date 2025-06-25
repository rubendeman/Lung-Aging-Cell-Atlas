setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

load('imm12_10.RData') #load integrated

DefaultAssay(immune.combined)<-'RNA'
immune.combined<-NormalizeData(immune.combined)

immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
immune.combined$Manuscript_Identity[immune.combined$Manuscript_Identity=='Baylor']<-'Dataset 1'
immune.combined$Manuscript_Identity[immune.combined$Manuscript_Identity=='IPF Cell Atlas']<-'Dataset 2'
immune.combined$agebin[immune.combined$agebin=='old']<-'Aged'
immune.combined$agebin[immune.combined$agebin=='young']<-'Young'

comp<-list(Epithelial=c('AT1','AT2','AT2B','AT2S','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage','DC'),Lymphoid=c('B','T','Mast','NK'))
cell.types<-unlist(comp)

immune.combined<-subset(immune.combined,subset=predicted.id %in% cell.types,invert=F)

comp_data<-immune.combined$predicted.id
for(i in names(comp)){comp_data[comp_data %in% comp[[i]]]<-i}
immune.combined<-AddMetaData(immune.combined,comp_data,col.name='comp')

#set.seed(1)
#fill_df <- immune.combined@meta.data %>% select('predicted.id','comp') %>% unique() %>% mutate(color = sample(scales::hue_pal()(25)))
#save(fill_df,file='fill_df.RData')
load('fill_df.RData')
fill_df$comp[fill_df$predicted.id=='DC']<-'Myeloid'

###########
# FIGURES #
###########
#UMAP
legends <- lapply(SplitObject(immune.combined, 'comp'), function(x) {
    comp<-x$comp[1]
  patchwork::wrap_elements(full = cowplot::get_legend(
    DimPlot(x,group.by='predicted.id',cols=fill_df$color[fill_df$comp==comp],order=rev(fill_df$predicted.id[fill_df$comp==comp]))+labs(color=comp)#+scale_fill_manual(name=comp,values=setNames(fill_df$color[fill_df$comp==comp],fill_df$predicted.id[fill_df$comp==comp]))
    ))
})
p1<-DimPlot(immune.combined,label=T,repel=T,group.by='predicted.id',cols=fill_df$color,order=rev(fill_df$predicted.id))+NoLegend()+
    theme(strip.background = element_rect(fill="red"), strip.text = element_text(size=15,face='bold',colour="white"),
    plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p1[[1]][["data"]][["tempvar"]]='Cell Type'
p0<-(legends[[3]]+legends[[2]])/(legends[[4]]+legends[[1]])/legends[[5]]
p2<-DimPlot(immune.combined,group.by='Manuscript_Identity',cols=c('#3E4E50','#FACFAD'),shuffle=T)+
    theme(legend.position=c(.75,.15), strip.background = element_rect(fill="green"), strip.text = element_text(size=15,face='bold',colour="white"),
    plot.title =element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p2[[1]][["data"]][["tempvar"]]='Dataset'
#p3<-FeaturePlot(immune.combined,features = 'age')+
p3<-DimPlot(immune.combined,group.by = 'agebin',cols=c('#91C4F2','#253C78'),shuffle=T)+
    theme(legend.position=c(0,.3), strip.background = element_rect(fill="blue"), strip.text = element_text(size=15,face='bold',colour="white"),
    plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p3[[1]][["data"]][["tempvar"]]='Age'
plot_grid(p1+facet_wrap(~tempvar), p0,p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar),ncol=2)

#MARKER GENES
DefaultAssay(immune.combined) = "RNA"
Idents(immune.combined)<-immune.combined$predicted.id
levels(immune.combined) <- c(rev(cell.types),levels(immune.combined)[!(levels(immune.combined) %in% cell.types)])
marker.genes<-c('AGER','SFTPC','KRT5','FOXJ1','CTSE','MUC5AC','PROX1','COL15A1','S100A3','IL7R','DKK2','ACKR1','PI16','PDGFRA','CDH11','ACTG2','TRPC6','S100A12','IL1B','ALOX5AP','GPR183','CD79A','IL32','TPSAB1','NKG7')
dot<-DotPlot(immune.combined, features = marker.genes, idents=cell.types)+ theme_minimal()+theme(axis.text.x = element_text(angle=90))+ RotatedAxis()+labs(y=NULL)

#DEMOGRAPHICS
tbl<-as.data.frame(table(immune.combined$orig.ident))
colnames(tbl)<-c('Subject','Age')
tbl$Age<-immune.combined$age[match(tbl$Subject,immune.combined$orig.ident)]
tbl<-tbl %>% mutate(Age=cut(Age,breaks=c(0,10,20,30,40,50,60,70,80,90,100)))
tbl$Age=substring(sub(',','-',tbl$Age),2,6)
tbl<-tbl %>% mutate(Sex=immune.combined$sex[match(tbl$Subject,immune.combined$orig.ident)])
tbl<-tbl %>% mutate(Source=immune.combined$Manuscript_Identity[match(tbl$Subject,immune.combined$orig.ident)])
tbl2=data.frame(Source=sample(immune.combined$Manuscript_Identity,5000))

#g1<-ggplot(tbl, aes(x="", y="", fill=Age)) + geom_bar(stat="identity", width=1) + coord_polar("y") + theme_void() + theme(legend.position=c(1.2,.5))
#g2<-ggplot(tbl, aes(x="", y="", fill=Sex)) + geom_bar(stat="identity", width=1) + coord_polar("y") + theme_void() + theme(legend.position=c(1.2,.5))
g12<-ggplot(tbl, aes(x=Age,fill=Sex))+ geom_histogram(stat="count")+theme_minimal()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_manual(values=c('#919098','#29E7CD'))+labs(y='Count',x='Age')
g3<-ggplot(tbl, aes(x="", y="", fill=Source)) + geom_bar(stat="identity", width=1) + coord_polar("y") + theme_void()+guides(fill=guide_legend(title="Source (Samples)"))+ theme(legend.position=c(1.2,.5))+scale_fill_manual(values=c('#3E4E50','#FACFAD'))
g4<-ggplot(tbl2, aes(x="", y="", fill=Source)) + geom_bar(stat="identity", width=1) + coord_polar("y") + theme_void()+guides(fill=guide_legend(title="Source (Cells)"))+ theme(legend.position=c(1.2,.5))+scale_fill_manual(values=c('#3E4E50','#FACFAD'))
plot_grid(dot,g12,plot_grid(g3,g4),ncol=1,rel_heights = c(2,1,1))+theme(plot.margin = unit(c(0, 1.5, 0, 0.5), "cm"))

#######################
### Cell Type Props ###
#######################
y1<-as.data.frame(cbind(immune.combined@meta.data$orig.ident,immune.combined@meta.data$predicted.id,immune.combined@meta.data$age))
y2=y1 %>% group_by(V1,V2) %>% summarise(V4=length(V1),V5=max(V3));
y3=y2 %>% group_by(V2) %>% summarise(V7=list(V1),V8=list(V4),V9=list(V5));
sids<-as.character(as.data.frame(table(immune.combined@meta.data$orig.ident))$Var1)
ind=apply(y3,1,function(x){unlist(x[3])[match(sids,unlist(x[2]))]}) #Create table of cell type versus sample
colnames(ind)<-y3$V2
rownames(ind)<-sids
demo_sub<-as.data.frame(lapply(comp,function(x){rowSums(ind[,x],na.rm=T)}))

#filtering
id_cond<-ind[,2]>5
idsub=rownames(ind)[id_cond]
ind=ind[id_cond,]
demo_sub=demo_sub[id_cond,]

demo=as.data.frame(cbind(count=apply(ind,1,function(x){sum(x,na.rm=T)}),age=immune.combined@meta.data$age[match(idsub,immune.combined@meta.data$orig.ident)]))

#ind<-ind/demo$count #Convert ind from raw counts to props
#or
for (i in 1:5){ind[,comp[[i]]]<-ind[,comp[[i]]]/demo_sub[,i]} #Convert ind from raw counts to sub props

p=list();
for (i in 1:5){
cormat2<-na.omit(as.data.frame(pivot_longer(as.data.frame(cbind(age=demo$age,ind[,comp[[i]]])),!age)))
vage<-cormat2$name
vage[cormat2$age<61]<-'Young'
vage[cormat2$age>60]<-'Aged'
cormat2$age<-vage
if(i==1){epidata=cormat2}
p[[i]]<-ggplot(cormat2, aes(y = value, x = name,fill=factor(age,levels=c('Young','Aged')))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Cell Type Proportion")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))+scale_fill_manual(values=c('#91C4F2','#253C78'))
}
plot_grid(plotlist=p)+theme(plot.margin = margin(1,1,1,1, "cm"))

wilcox.test(epidata$value[epidata$name=='AT2'&epidata$age=='Aged'],epidata$value[epidata$name=='AT2'&epidata$age=='Young'],alternative='less')

###################
### Aging Genes ###
###################
#PRINT TO CONSOLE
load('humandiff_1_12_24.RData')
res<-datalistlist[[1]]$gCap
noquote(head(rownames(res),200))
upgenes<-rownames(res)[res$avg_log2FC>0][1:20]
dngenes<-rownames(res)[res$avg_log2FC<0][1:20]

#VOLCANO
res0<-datalistlist[[1]]$gCap
res<-res0
res<-res[(!grepl("^MT", rownames(res))),]
#res<-res[rownames(res)%in%gnames,]
with(res,plot(avg_log2FC, -log10(p_val_adj), pch = 19, xlab='log2FC', ylab=expression('-log'[10]*'(p)')))
abline(v=0)
abline(h=-log10(0.05))
res2=res[res$p_val_adj<0.05,]
xind=abs(res2$avg_log2FC)>1#|res2$p_val_adj<0.000000000000000000000000001
text(res2$avg_log2FC[xind],-log10(res2$p_val_adj[xind]),rownames(res2)[xind],cex=0.5,pos=1)

#GENE LISTS
feat_check=list()
l=load('all_sen_lists.RData')
#feat_check[['sen']]=unlist(c(senmayolist,cell,fridman,segura,purcell))
feat_check[['igf']] <- c('PAPPA','IGF1','IGF1R','IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7')
feat_check[['focal']] <- c('ITGB1','FLNA','ITGA1','ITGA2','PARVA','CACNB2') #Focal adhesion in fibroblasts
feat_check[['at2']] <- c('SFTPA1','SFTPA2','SFTPB','SFTPC','ABCA3','HHIP','LAMP3') #Surfactants in AT2
feat_check[['ifn']] <- c('STAT1','STAT2','IFIT1','JAK2','NLRC5','TXK','PARP9') #Interferon genes in AM
feat_check[['mtdna']] <- c('MT-ATP8','MT-ATP6','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND5','MT-ND6') #mtDNA
feat_check[['cap']]=c('VIPR1','FENDRR','FOXF1','CA4','F2RL3','IL7R','APLNR','GPIHBP1','FCN3','EDN1','CA4','APLN','EDNRB','HPGD','DKK2','GJA5','SERPINE2','HEY1','EFNB2','ACKR1','PRSS23','VWF') #LGEA
feat_check[['ubiq']]=c('ABLIM1','MIB1','NPEPPS','UBXN7')
feat_check[['mapk']]=c('ANXA4','MAPK8','MAP3K20','RABGAP1L','RAPGEF2')

#VIOLIN PLOTS
Idents(immune.combined)<-immune.combined$predicted.id
sub.obj=subset(immune.combined,subset=predicted.id %in% c('gCap','Aerocyte'))
sub.obj$agebin<-factor(sub.obj$agebin,levels=c('Young','Aged'))

VlnPlot(sub.obj,features=feat_check[['mtdna']],split.by='agebin')

VlnPlot(sub.obj,features=c('VIPR1','FOXF1','CA4','IL7R','FCN3','HPGD'),split.by='agebin')
VlnPlot(sub.obj,features=c('VIPR1','IL7R','FCN3','HPGD'),split.by='agebin',cols=c('#253C78','#91C4F2'),ncol=4)&theme(axis.title.x = element_blank())&stat_compare_means(label='p.signif')

#DOTPLOT
Idents(immune.combined)<-immune.combined$predicted.id
a<-DotPlot(immune.combined,idents=cell.types,features=feats,group.by='predicted.id',split.by='agebin',cols='Spectral',scale=F)
a$data$id<-sub("_", " ",a$data$id)
a

#BOXPLOT SAMPLE MEANED EXPRESSION
cell.types='AT2'

genes = feats = 'SFTPC'
cellorgene<-'gene'
DefaultAssay(immune.combined)<-'RNA'
datalist = list(); #Empty list for results
type=cell.types
#Iterate through cell types
for (i in 1:length(type)){
celltype=type[i];
Idents(immune.combined)=immune.combined@meta.data$predicted.id
tmp=subset(immune.combined, subset=(predicted.id %in% celltype));
#DotPlot
Idents(tmp)<-tmp@meta.data$orig.ident
a<-DotPlot(tmp,features=genes,scale=F)
mtc<-match(a$data$id,tmp@meta.data$orig.ident)
a$data$id<-tmp@meta.data$agebin[mtc]
newdata<-a$data #[a$data$avg.exp!=0,]
newdata2=newdata %>% group_by(id) %>% summarise(avg=mean(avg.exp),pct=mean(pct.exp))

#or AverageExpression
#avgexp<-as.data.frame(AverageExpression(tmp,assays='RNA',genes,group.by = c('orig.ident','agebin'))$RNA)
#avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))})) #IMPORTANT; SCALE to 0-1
#avgexp=data.frame(gene=genes,avgexp)
#newdata=pivot_longer(avgexp,!gene)
#newdata=as.data.frame(as.matrix(within(newdata, name<-data.frame(do.call('rbind', strsplit(as.character(name), '_', fixed=TRUE))))))
#colnames(newdata)<-c('features.plot','sample','id','avg.exp')
#newdata$avg.exp=as.numeric(newdata$avg.exp)

datalist[[i]]=newdata
}
names(datalist)=type

datalist2<-bind_rows(datalist, .id = "column_label")
datalist2<-datalist2[datalist2$avg.exp!=0,] #IMPORTANT; DROP 0 EXPRESSION
if(cellorgene=='gene'){datalist2$column_label=datalist2$features.plot}
featbox<-ggplot(datalist2, aes(x=column_label, y=avg.exp,fill=factor(id,levels=c('Young','Aged')))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Gene Expression")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))
featbox

#BIG MAP BY SAMPLE (NEW)
cell.types=c('Aerocyte','gCap')
ht=NULL
for (i in cell.types){
    feats=c(upgenes,dngenes)
    rowsplit=c(rep('Increased',length(upgenes)),rep('Decreased',length(dngenes)))
    noquote(feats)
    #Downsample
    cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i))#&Manuscript_Identity=='Dataset 1')) #BAYLOR ONLY
    cbmc.small<-NormalizeData(cbmc.small)
    cbmc.small<-ScaleData(cbmc.small,features=feats,split.by = 'Manuscript_Identity')
    Idents(cbmc.small)<-paste(cbmc.small$predicted.id,cbmc.small$agebin)
    avgexp<-AverageExpression(cbmc.small,assays='RNA',layer='scale.data',feats,group.by = c('predicted.id','orig.ident','agebin'))
    avgexp<-as.matrix(avgexp$RNA)
    
    #Scale Matrix
    #avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))}))
    avgexp[is.na(avgexp)]<-0
    
    #Create Annotations
    ann <- data.frame(i,sub('_.*','',colnames(avgexp)),sub('.*_','',colnames(avgexp)))
    colnames(ann) <- c('Cell Type','Sample','Age')
    ann<-ann %>% mutate(Dataset=cbmc.small$Manuscript_Identity[match(ann$Sample,cbmc.small$orig.ident)])
    ann<-ann[,c(1,3,4,2)]
    #ann <- data.frame(sub('_.*','',colnames(avgexp)),cut(cbmc.small$age[match(gsub('.*_(.+)_.*','\\1',colnames(avgexp)),cbmc.small$orig.ident)],c(0,40,60,100),c('young','mid','old')))
    ord=order(ann$Dataset)
    ann=ann[ord,]
    avgexp=avgexp[,ord]

    hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'))
    names(hcols$'Cell Type')=i
    colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))
    rowfeats=NULL #unlist(feat_check)
    ha = rowAnnotation(foo = anno_mark(at = which(rownames(avgexp) %in% rowfeats), labels = rownames(avgexp)[rownames(avgexp)%in%rowfeats]))

    col_fun = circlize::colorRamp2(c(-1,0, 1), c("purple", "black", "#FFFF00"))
    ht0=Heatmap(na.omit(avgexp), name = "Expression",  column_split=data.frame(factor(ann$'Cell Type'),factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
        cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
        show_row_dend = F, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
        show_column_names = F, show_row_names=T, row_split=factor(rowsplit,levels=c('Increased','Decreased')), use_raster = F)
    ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 8, 2), "mm"),merge_legend=T,ht_gap=unit(-1,'cm'))

#BIG MAP BY SAMPLE (OLD VERSION)
cell.types=c('Aerocyte','gCap')
ht=NULL
for (i in cell.types){
feats=c(upgenes,dngenes)
rowsplit=c(rep('Increased',length(upgenes)),rep('Decreased',length(dngenes)))
noquote(feats)
#Downsample
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i))#&Manuscript_Identity=='Dataset 1')) #BAYLOR ONLY
cbmc.small<-NormalizeData(cbmc.small)
Idents(cbmc.small)<-paste(cbmc.small$predicted.id,cbmc.small$agebin)
avgexp<-AverageExpression(cbmc.small,assays='RNA',layer='data',feats,group.by = c('predicted.id','orig.ident','agebin'))
avgexp<-as.matrix(avgexp$RNA)

#Scale Matrix
#avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))}))
#or scale datasets independently. 7/10/24
allind=as.data.frame(table(paste(immune.combined$orig.ident,immune.combined$Manuscript_Identity)))
allind=data.frame(do.call('rbind', strsplit(as.character(allind$Var1),' ',fixed=TRUE)))
ind1=na.omit(match(allind$X1[allind$X3==1],stringr::str_split_fixed(colnames(avgexp),"_",2)[,1]))
ind2=na.omit(match(allind$X1[allind$X3==2],stringr::str_split_fixed(colnames(avgexp),"_",2)[,1]))
#avgexp[,ind1]=t(apply(avgexp[,ind1],1,function(x){(x-min(x))/(max(x)-min(x))}))
#avgexp[,ind2]=t(apply(avgexp[,ind2],1,function(x){(x-min(x))/(max(x)-min(x))}))
avgexp=cbind(t(apply(avgexp[,ind1],1,function(x){(x-min(x))/(max(x)-min(x))})),t(apply(avgexp[,ind2],1,function(x){(x-min(x))/(max(x)-min(x))})))
avgexp[is.na(avgexp)]<-0

#Create Annotations
if(length(i)==1){
    ann <- data.frame(i,sub('_.*','',colnames(avgexp)),sub('.*_','',colnames(avgexp)))
    colnames(ann) <- c('Cell Type','Sample','Age')
    ann<-ann %>% mutate(Dataset=cbmc.small$Manuscript_Identity[match(ann$Sample,cbmc.small$orig.ident)])
    ann<-ann[,c(1,3,4,2)]
#ann <- data.frame(sub('_.*','',colnames(avgexp)),cut(cbmc.small$age[match(gsub('.*_(.+)_.*','\\1',colnames(avgexp)),cbmc.small$orig.ident)],c(0,40,60,100),c('young','mid','old')))
}else{
}
hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))
rowfeats=NULL #unlist(feat_check)
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avgexp) %in% rowfeats), labels = rownames(avgexp)[rownames(avgexp)%in%rowfeats]))

col_fun = circlize::colorRamp2(c(0.2, 1), c("black", "#FFFF00"))
#col_fun = circlize::colorRamp2(c(0,0.5, 1), c("purple", "black", "#FFFF00"))
ht0=Heatmap(na.omit(avgexp), name = "Expression",  column_split=data.frame(factor(ann$'Cell Type'),factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
    cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
    show_row_dend = F, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
    show_column_names = F, show_row_names=T, row_split=factor(rowsplit,levels=c('Increased','Decreased')), use_raster = F)
    ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 8, 2), "mm"),merge_legend=T,ht_gap=unit(-1,'cm'))

#BIG MAP BY CELL
ht=NULL
rowsplit=NULL
cell.types=c('Aerocyte','gCap')
for (i in cell.types){
feats=c(upgenes,dngenes)
rowsplit=c(rep('Increased',length(upgenes)),rep('Decreased',length(dngenes)))
noquote(feats)
#Downsample
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% i)) #&Manuscript_Identity=='Dataset 1'))
Idents(cbmc.small)<-cbmc.small$orig.ident
table(cbmc.small$orig.ident)
cbmc.small <- subset(cbmc.small, downsample = 100)
Idents(cbmc.small)<-'avgexp'
cbmc.small <- subset(cbmc.small, downsample = 5000)
avgexp=cbmc.small@assays$RNA@data[feats,] %>% as.matrix()

#Scale Matrix
avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))}))
avgexp[is.na(avgexp)]<-0

#Create Annotations
ann <- data.frame(cbmc.small$predicted.id, cbmc.small$agebin, cbmc.small$Manuscript_Identity, cbmc.small$orig.ident)
colnames(ann) <- c('Cell Type','Age','Dataset','Sample')
hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=c(Aged='#91C4F2',Young='#253C78'),Dataset=c('Dataset 1'='#3E4E50','Dataset 2'='#FACFAD'))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(T,T,T,F))
rowfeats=unlist(feat_check)
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avgexp) %in% rowfeats), labels = rownames(avgexp)[rownames(avgexp)%in%rowfeats]))

#col_fun = circlize::colorRamp2(c(-1, 0, 1), c("purple", "black", "#FFFF00"))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht0=Heatmap(avgexp, name = "Expression",  column_split=data.frame(factor(ann$'Cell Type'),factor(ann$'Age',levels=c('Young','Aged'))),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = T,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=F, row_split=rowsplit, use_raster = F)
            ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 8, 2), "mm"),merge_legend=T)

##################
### Senescence ###
##################
range(immune.combined@meta.data$SenMayo1)
#Density plot to determine threshold (i.e. High SenMayo cells are above this threshold)
plot(density(immune.combined@meta.data$SenMayo1),main='SenMayo Score Distribution',frame.plot=F,xlim=c(-0.2,0.6),xaxt="n",ylim=c(0,7))
thresh=0.25
abline(v=thresh,col="red")
axis(1,pos=0)

sum(immune.combined@meta.data$SenMayo1>thresh)
ages3  <- vector(mode='character',length=length(immune.combined@meta.data$SenMayo1))
for (x in 1:length(immune.combined@meta.data$SenMayo1)) {
    val=immune.combined@meta.data$SenMayo1[x]; if (val<thresh){ages3[x]="Low SenMayo";} else {ages3[x]="High SenMayo";}
}
table(ages3)
immune.combined <- AddMetaData(object = immune.combined, metadata = ages3, col.name = 'senbin')

#Senescence Cell Type Scores
data=data.frame(ct=immune.combined$predicted.id,sm=immune.combined$SenMayo1)
data$ct=factor(data$ct,levels=cell.types)
ggplot(data, aes(x = ct,y = sm)) + geom_boxplot() + coord_flip()+theme_minimal()+xlab('Cell Type')+ylab('SenMayo Score')

#Senescence VS age
genes = feats = 'SenMayo1'
cellorgene<-'cell'
DefaultAssay(immune.combined)<-'RNA'
datalist = list(); #Empty list for results
type=cell.types
#Iterate through cell types
for (i in 1:length(type)){
celltype=type[i];
Idents(immune.combined)=immune.combined@meta.data$predicted.id
tmp=subset(immune.combined, subset=(predicted.id %in% celltype));
#DotPlot
Idents(tmp)<-tmp@meta.data$orig.ident
a<-DotPlot(tmp,features=genes,scale=F)
mtc<-match(a$data$id,tmp@meta.data$orig.ident)
a$data$id<-tmp@meta.data$agebin[mtc]
newdata<-a$data #[a$data$avg.exp!=0,]
newdata2=newdata %>% group_by(id) %>% summarise(avg=mean(avg.exp),pct=mean(pct.exp))
datalist[[i]]=newdata
}
names(datalist)=type

datalist2<-bind_rows(datalist, .id = "column_label")
datord=as.data.frame(datalist2 %>% group_by(column_label) %>% summarise(avg=mean(avg.exp)))
datord=datord$column_label[order(datord$avg,decreasing = T)] #USE FOR LEVELS ARG X AXIS TO ORDER BY MAGNITUDE
if(cellorgene=='gene'){datalist2$column_label=datalist2$features.plot}
featbox<-ggplot(datalist2, aes(x=factor(column_label,levels=cell.types), y=avg.exp,group=interaction(factor(id,levels=c('Young','Aged')),column_label),fill=factor(column_label,levels=datord))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "SenMayo Score")+geom_point(position=position_jitterdodge(jitter.width=0.1))+
    scale_fill_manual(values=fill_df$color[match(datord,fill_df$predicted.id)])+NoLegend()
featbox

datalist2<-datalist2 %>% mutate(color=fill_df$color[match(datalist2$column_label,fill_df$predicted.id)])
d3<-datalist2 %>%  mutate(class=factor(immune.combined$comp[match(column_label,immune.combined$predicted.id)],levels=names(comp))) %>% group_by(class)
d3$column_label=factor(d3$column_label,levels=cell.types)
featbox<-ggplot(d3, aes(x=column_label, y=avg.exp,group=interaction(factor(id,levels=c('Young','Aged')),column_label))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "SenMayo Score")+geom_point(aes(color=column_label),position=position_jitterdodge(jitter.width=0.1))+facet_wrap(~ class,nrow=1,scales='free_x')+
    NoLegend()+scale_color_manual(values=unique(d3$color))+stat_compare_means(label='p.signif')

#SENESCENCE MAIN FIGURE
p1<-DimPlot(immune.combined,raster=F,pt.size=0.01,label=T,repel=T,group.by='predicted.id',cols=fill_df$color,order=rev(fill_df$predicted.id))+NoLegend()+labs(title=NULL,x="UMAP 1", y = "UMAP 2")
p2<-DimPlot(immune.combined,raster=F,pt.size=0.01,group.by='senbin',order='High Senescence Score',cols = c("#4A306D", "#FAB2EA"))+NoLegend()+labs(title=NULL, x="UMAP 1", y = "UMAP 2")
p3<-VlnPlot(immune.combined, features = c("CDKN1A"), group.by = "senbin",  pt.size = 0, cols=c('#FAB2EA','#4A306D'))+theme(line=element_blank(), axis.text.x =element_blank(), axis.title = element_blank())+ggtitle('CDKN1A')

plot_grid(p1, p2,ncol=2)

rightup=plot_grid(NULL,p3,NULL,ncol=3,rel_widths = c(0.5,2,0.5))
plot_grid(rightup,featbox,ncol=1,rel_heights = c(1,2))+ theme(plot.margin = unit(c(0, 0, 0.5, 1), "cm"))

#SENESCENCE MODULES UMAP PLOTS
gene_names=c('2','15','19','31')
gg_Fig <- FeaturePlot(immune.combined, features = c("X21","X151","X191","X311"),ncol=2,raster=F)
gg_Fig <- lapply((1:4), function(x) { gg_Fig[[x]] + labs(title=gene_names[x])+labs(x="UMAP 1", y = "UMAP 2")})
CombinePlots(gg_Fig)

#BIG MAP SENMAYO MODULES
load('geneInfo0.RData')
sencell=data.frame(Gene=unlist(read.table('listsenmayo.txt'))) #mapping high senescence score cells
mods=c(2,15,19,20,28,31)
sencell=sencell %>% mutate(Color=geneInfo0$moduleColor[match(sencell$Gene,geneInfo0$geneSymbol)])
sencell=sencell[sencell$Color %in% mods,]

Idents(immune.combined)<-immune.combined$predicted.id
downsample=subset(immune.combined,downsample=100)
avgexp=downsample@assays$RNA@data[sencell$Gene, ] %>% as.matrix()
avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht=Heatmap(na.omit(avgexp), name = "Expression", column_split = factor(downsample$predicted.id,levels=cell.types), cluster_columns = F, show_column_dend = FALSE,
            cluster_column_slices = F, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), row_split = factor(sencell$Color,levels=mods),cluster_rows = T,
            cluster_row_slices=F, show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 6), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, use_raster = TRUE, raster_quality = 4)
draw(ht, padding = unit(c(2, 2, 24, 10), "mm"))

###############################
# Senescence-correlated genes #
###############################
#read glmmTMB output
cell.types<-c('AT2','gCap')
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/senage/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})
names(datal)<-cell.types

#venn
library(VennDiagram)
overrideTriple=T
datal2=lapply(datal,function(x){x$Gene[x$p_val_adj<0.05&x$Beta1>0]})
v1=venn.diagram(x=datal2,filename=NULL,fill=fill_df$color[match(cell.types,fill_df$predicted.id)])
plot_grid(v1)

#find unique genes
library(gplots)
ItemsList <- venn(datal2, show.plot=FALSE)
attributes <- attr(ItemsList,"intersections")
data_unique<-lapply(cell.types,function(x){datal[[x]][datal[[x]]$Gene %in% attributes[[x]],]})

#venn anti-correlated
datal2=lapply(datal,function(x){x$Gene[x$p_val_adj<0.05&x$Beta1<0]})
v2=venn.diagram(x=datal2,filename=NULL,fill=fill_df$color[match(cell.types,fill_df$predicted.id)])
plot_grid(v2)
plot_grid(v1,v2,ncol=1)

#modules
load('geneInfo0.RData')
which(geneInfo0$geneSymbol[geneInfo0$moduleColor %in% c(19,20)] %in% datal2[['SMC']])

#hmap
cell.types<-c('AT2','gCap')
datal=data_unique
names(datal)<-cell.types
datal=lapply(datal,function(x){x[order(x$Beta1,decreasing=T)[1:10],]})
sen2=bind_rows(datal, .id = "column_label")
ht=NULL
feats=sen2$Gene
rowsplit=sen2$column_label
cbmc.small<-subset(immune.combined,subset=(predicted.id %in% cell.types))
Idents(cbmc.small)<-cbmc.small$predicted.id
cbmc.small <- subset(cbmc.small, downsample = 100)
#CREATE MATRIX
order=order(cbmc.small$SenMayo1,decreasing=F)
avgexp=cbmc.small@assays$RNA@data[feats,order] %>% as.matrix()
#CREATE ANNOT
ann <- data.frame(cbmc.small$predicted.id[order], cbmc.small$orig.ident[order],cbmc.small$SenMayo1[order])
colnames(ann) <- c('Cell Type','Sample','SenScore')
col_fun = circlize::colorRamp2(c(0, 25), c("grey", "blue"))
hcols=list('Cell Type'=fill_df$color[match(cell.types,fill_df$predicted.id)],Beta1=col_fun)
names(hcols$'Cell Type')=cell.types
colAnn <- HeatmapAnnotation(df = ann,col=hcols,which = 'col',gap = unit(1, 'mm'),show_legend=c(F,F,T))
ann2 <- data.frame(ct=rowsplit,Beta1=sen2$Beta1)
colnames(ann2) <- c('Cell Type','Beta1')
rowfeats=NULL
ha = rowAnnotation(df=ann2,col=hcols,foo = anno_mark(at = which(rownames(avgexp) %in% rowfeats),labels = rownames(avgexp)[rownames(avgexp)%in%rowfeats]),annotation_legend_param=list(Beta1=list(at=seq(0,25,5))))

avgexp=t(apply(avgexp,1,function(x){(x-min(x))/(max(x)-min(x))}))
#col_fun = circlize::colorRamp2(c(-1, 0, 1), c("purple", "black", "#FFFF00"))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht=Heatmap(avgexp, name = "Expression",  column_split=data.frame(factor(ann$'Cell Type')),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F, column_title=NULL, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = F,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=T, row_split=rowsplit, use_raster = F)
draw(ht, padding = unit(c(8, 2, 2, 8), "mm"),merge_legends=T)

##################
# P16 Expression #
##################
datalist = list();
type=cell.types
#Iterate through cell types
for (i in 1:length(type)){
celltype=type[i];
Idents(immune.combined)=immune.combined@meta.data$predicted.id
tmp=subset(immune.combined, idents = celltype);
Idents(tmp)=tmp@meta.data$agebin;
DefaultAssay(tmp)<-'RNA'
agegenes <- FindMarkers(tmp, features='CDKN2A', ident.1 = "old", ident.2 = "young", logfc.threshold = 0, min.pct = 0,verbose = FALSE);
datalist[[i]]=na.omit(agegenes);
}
names(datalist)=type

#P16 Prop Cor
DefaultAssay(immune.combined)<-'RNA'
Idents(immune.combined)<-immune.combined@meta.data$orig.ident
a<-DotPlot(immune.combined,features='CDKN2A')
a$data$id<-immune.combined@meta.data$age[match(a$data$id,immune.combined@meta.data$orig.ident)]
View(a$data)
pctexp<-a$data$pct.exp
aage<-a$data$id
cor(pctexp,aage)

data=data.frame(Age=aage,Score=pctexp)
data$Age[data$Age>60]<-'61-80'
data$Age[data$Age<61&data$Age>40]<-'41-60'
data$Age[data$Age<41&data$Age>20]<-'21-40'
data$Age[data$Age<21]<-'0-20'
p<-ggplot(data, aes(x=Age, y=Score)) + geom_boxplot() + labs(x='Age',y='% Cells Expressing p16')
p + theme_classic()+theme(text = element_text(size=15)) + scale_y_continuous(limits=c(0,10))

#############################
### TARGETED GENE AGE COR ###
#############################
#Corrected by #cells each subject. Later, let's also look at # sen cells per subject VS age
sencell=c('CDKN1A','HHIP','GLB1','NEU1','IGFBP7','QDPR');
corv = list()
pv=list()
nn=list()
for (i in sencell){
rnadat=immune.combined@assays$RNA@data[i,] %>% as.matrix()
m=data.frame(immune.combined@meta.data$orig.ident,immune.combined@meta.data$predicted.id,as.numeric(immune.combined@meta.data$age),rnadat)
thresh1=0.1
thresh2=0.6
o=m %>% group_by(immune.combined.meta.data.predicted.id,immune.combined.meta.data.orig.ident)  %>% summarise(V1=mean(as.numeric.immune.combined.meta.data.age.),V2=mean(rnadat),V3=length(rnadat))
p=o %>% group_by(immune.combined.meta.data.predicted.id) %>% summarise(V4=cor(V1,V2,method="spearman"),V5=corPvalueStudent(cor(V1,V2,method="spearman"),sum(V3)),V6=sum(V2>thresh1)/length(V2))
corv[[i]]=p[,2]
pv[[i]]=p[,3]
nn[[i]]=p[,4]
}
corv=do.call("cbind",corv)
rownames(corv)=p$immune.combined.meta.data.predicted.id
colnames(corv)=sencell
pv=do.call("cbind",pv)
rownames(pv)=p$immune.combined.meta.data.predicted.id
colnames(pv)=sencell
nn=do.call("cbind",nn)
nnn=nn<thresh2
corv2=corv
corv2[nnn]=NA
#corv2[is.na(corv2)]=-1
col_fun = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#FF00FF","black" ,"#FFFF00"))
#corv2[corv2==-1]=NA
h=Heatmap((corv2[c(-3,-4),]), col=col_fun, name = "Age Correlation", show_column_dend = FALSE, column_title_rot = 45, cluster_columns = FALSE, 
cluster_rows = FALSE,show_row_dend = FALSE, row_names_gp = gpar(fontsize = 12),row_names_side = "left")
draw(h, padding = unit(c(10, 10, 2, 10), "mm"))

immune.combined$predicted.id=factor(immune.combined$predicted.id,levels=cell.types)
Idents(immune.combined)<-immune.combined$predicted.id
VlnPlot(immune.combined,features=c('HHIP','GLB1','CDKN1A','IGFBP7'),pt.size=0,ncol=2)&theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, hjust=1))

