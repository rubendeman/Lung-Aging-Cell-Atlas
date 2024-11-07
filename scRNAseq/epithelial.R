setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ComplexHeatmap)

#EPITHELIUM: USE THE CURRENT SEURAT OBJECT... OR...
Idents(immune.combined)=immune.combined@meta.data$predicted.id
epi <- subset(immune.combined, idents=comp[['Epithelial']])

DefaultAssay(epi)<-"integrated"
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:20)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:20)
epi <- FindClusters(epi, resolution = 1)
DimPlot(epi,label=T,raster=F)+DimPlot(epi,group.by='predicted.id',raster=F,label=T)

#...LOAD the saved one
load('epi_4_11_24.RData')
DimPlot(epi)
epi$predicted.id[epi$predicted.id=='AT2']='AT2B'='SPCʰⁱᵍʰ AT2'
epi$predicted.id[epi$integrated_snn_res.0.5==9]='AT2S'='SPCˡᵒʷ AT2'
epi$agebin[epi$agebin=='old']='Aged'
epi$agebin[epi$agebin=='young']='Young'

#find at2s markers
#epi=subset(epi,subset=Manuscript_Identity=='Dataset 1')
Idents(epi)<-epi$predicted.id
DefaultAssay(epi)<-'RNA'
mrk=FindMarkers(epi,ident.1='SPCˡᵒʷ AT2',ident.2='SPCʰⁱᵍʰ AT2')
noquote(rownames(mrk)[mrk$avg_log2FC>0][1:500])

genes=c('SFTPC','SFTPD','HHIP','MKI67','HMMR','TOP2A','ERBB4','TNIK','TCF12','FOXP1','STAT3','YAP1','TEAD1','AXIN1','APC','CTNNB1','LRP5','WNT10A')
mrk[genes,]
VlnPlot(epi,c('SFTPA1','SFTPA2','SFTPB','SFTPC'),pt.size=0,idents=c('AT1','SPCʰⁱᵍʰ AT2','SPCˡᵒʷ AT2'),cols=c('#ED8141','#5BB300','red'),ncol=4)&scale_y_continuous(limits=c(0,NA))&theme(axis.title.x = element_blank())
VlnPlot(epi,c('CTNNB1','GSK3B','TCF12','ERBB4','JAK1','STAT3','YAP1','STK3'),pt.size=0,idents=c('AT1','SPCʰⁱᵍʰ AT2','SPCˡᵒʷ AT2'),cols=c('#ED8141','#5BB300','red'),ncol=4)&scale_y_continuous(limits=c(0,NA))&theme(axis.title.x = element_blank())

library(glmmTMB)
soup.subset=subset(epi,subset=predicted.id %in% c('SPCˡᵒʷ AT2','SPCʰⁱᵍʰ AT2'))
which(genes %in% rownames(soup.subset@assays$RNA@counts))
subject <- soup.subset@meta.data$orig.ident
sf <- soup.subset@meta.data$nCount_RNA
age <- as.factor(soup.subset$predicted.id)
count_mat <- soup.subset@assays$RNA@counts[genes, ]
res.t <- apply(count_mat, 1, function(gene){
tmp.df <- data.frame(y=gene, norm_sf=sf, age=as.numeric(age), sub=as.numeric(factor(subject)))
  f2 <- try(glmmTMB(y ~ age + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = nbinom2))
  res2 <-try(c(
    summary(f2)$coefficients$cond[,"Estimate"],
    summary(f2)$coefficients$cond["age","Pr(>|z|)"]))
  if('try-error' %in% class(res2)){
    res2 <- rep(NA,3)
  }
  names(res2) <- c("Beta0","Beta1","P-Beta")
  return(res2)
})

#new plots
load('fill_df.RData')
epi_df=fill_df[fill_df$comp=='Epithelial',]
epi_df[4,1]='SPCʰⁱᵍʰ AT2'
epi_df=rbind(epi_df,c('SPCˡᵒʷ AT2','Epithelial','red'))
p1<-DimPlot(epi,label=T,repel=T,group.by='predicted.id',cols=epi_df$color,order=rev(epi_df$predicted.id),shuffle=T)+NoLegend()+
        theme(strip.background = element_rect(fill="red"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p1[[1]][["data"]][["tempvar"]]='Cell Type'
p2<-DimPlot(epi,group.by='Manuscript_Identity',shuffle=T,cols=c('#3E4E50','#FACFAD'))+
        theme(legend.position=c(.05,.85),
            strip.background = element_rect(fill="green"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title =element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p2[[1]][["data"]][["tempvar"]]='Dataset'
p3<-DimPlot(epi,group.by = 'agebin',cols=c('#91C4F2','#253C78'),shuffle=T)+
        theme(legend.position=c(.05,.8),
            strip.background = element_rect(fill="blue"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p3[[1]][["data"]][["tempvar"]]='Age'
p4<-FeaturePlot(epi,features = 'SFTPC')+
        theme(legend.position=c(.05,.8),
            strip.background = element_rect(fill="purple"),
              strip.text = element_text(size=15,face='bold',colour="white"),
              plot.title=element_blank(), line=element_blank(), axis.text =element_blank(), axis.title = element_blank())
p4[[1]][["data"]][["tempvar"]]='SFTPC'
#plot_grid(p1+facet_wrap(~tempvar), p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar), p4+facet_wrap(~tempvar),ncol=2)
plot_grid(p1+facet_wrap(~tempvar), p2+facet_wrap(~tempvar), p3+facet_wrap(~tempvar),ncol=3,byrow=T)

#HHIP Plots
Idents(epi)<-epi$predicted.id
epi$predicted.id<-factor(epi$predicted.id,levels=epi_df$predicted.id)
p5<-VlnPlot(epi, group.by='predicted.id',features = c("HHIP"),cols=epi_df$color,pt.size=0)+ggtitle('HHIP')+theme(legend.position = 'none')
p6<-FeaturePlot(epi,features='HHIP',label=T,repel=T,raster = F) + labs(x="UMAP 1", y = "UMAP 2")
p7<-FeaturePlot(epi,features='SFTPC',label=T,repel=T,raster = F) + labs(x="UMAP 1", y = "UMAP 2")
plot_grid(p7,p6,p5,ncol=3,byrow=T)

#First, let's plot props
y1<-as.data.frame(cbind(epi@meta.data$orig.ident,as.character(epi@meta.data$predicted.id),epi@meta.data$age))
y2=y1 %>% group_by(V1,V2) %>% summarise(V4=length(V1),V5=max(V3));
y3=y2 %>% group_by(V2) %>% summarise(V7=list(V1),V8=list(V4),V9=list(V5));
sids<-as.character(as.data.frame(table(epi@meta.data$orig.ident))$Var1)
ind=apply(y3,1,function(x){unlist(x[3])[match(sids,unlist(x[2]))]}) #Create table of cell type versus sample
colnames(ind)<-y3$V2
rownames(ind)<-sids
ind[is.na(ind)]=0

#First, let's filter out samples with very few AT2
id_cond<-ind[,1]>5#&ind[,2]>5
idsub=rownames(ind)[id_cond]
epi_sub=subset(epi,subset=(orig.ident %in% idsub)&(predicted.id %in% c('SPCˡᵒʷ AT2','SPCʰⁱᵍʰ AT2')))
ind=ind[id_cond,]

#Now, cell type props
demo=as.data.frame(cbind(count=rowSums(ind),age=epi@meta.data$age[match(rownames(ind),epi@meta.data$orig.ident)]))
indprop<-ind/demo$count
cormat2<-na.omit(as.data.frame(pivot_longer(as.data.frame(cbind(age=demo$age,indprop)),!age)))
cormat2$age<-cut(cormat2$age,c(0,60,100),c('Young','Aged'))
h0=ggplot(cormat2, aes(y = value, x = name,fill=factor(age,levels=c('Young','Aged')))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Cell Type Proportion")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))+scale_fill_manual(values=c('#91C4F2','#253C78'))

#Plot %AT2S young versus old
newdata=data.frame(epi_sub@meta.data$predicted.id, age=as.numeric(epi_sub@meta.data$age), sample=epi_sub@meta.data$orig.ident);
newdata=newdata %>% group_by(sample) %>% summarise(V1=sum(epi_sub.meta.data.predicted.id=='SPCˡᵒʷ AT2',na.rm=TRUE)/length(epi_sub.meta.data.predicted.id),age=median(age))
newdata$age=cut(newdata$age,c(0,60,90),c('Young','Aged'))
newdata=newdata[newdata$V1>0,]
h1<-ggplot(newdata,aes(x=age, y=V1)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
labs(x='Age',y='% SPCˡᵒʷ AT2')+theme_classic()+theme(text = element_text(size=15)) +scale_y_continuous(limits=c(0,1))
t.test(newdata$V1[newdata$age=='Aged'],newdata$V1[newdata$age=='Young'])

#Plot % expressing HHIP in AT2B vs AT2S
newdata=data.frame(cell=epi_sub@meta.data$predicted.id, sample=epi_sub@meta.data$orig.ident, gene=epi_sub@assays$RNA@data['HHIP',]>0)
newdata=newdata %>% group_by(sample,cell) %>% summarise(V1=sum(gene)/length(gene))
h2<-ggplot(newdata,aes(x=cell, y=V1)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
labs(x='Cell Type',y='% Expressing HHIP')+theme_classic()+theme(text = element_text(size=15)) +scale_y_continuous(limits=c(0,1))
plot_grid(h1,h2,ncol=2)

#PLOT OF SURFACTANTS
genes = feats = c('SFTPA1','SFTPA2','SFTPB','SFTPC')
cellorgene<-'gene'
DefaultAssay(epi_sub)<-'RNA'
datalist = list(); #Empty list for results
type='AT2'
#Iterate through cell types
for (i in 1:length(type)){
celltype=type[i];
Idents(epi_sub)=epi_sub@meta.data$predicted.id
tmp=subset(epi_sub, subset=(predicted.id %in% c('SPCˡᵒʷ AT2','SPCʰⁱᵍʰ AT2')));
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
datalist2$id[datalist2$id=='young']='Young'
datalist2$id[datalist2$id=='old']='Aged'
if(cellorgene=='gene'){datalist2$column_label=datalist2$features.plot}
featbox<-ggplot(datalist2, aes(x=column_label, y=avg.exp,fill=factor(id,levels=c('Young','Aged')))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Gene Expression")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))+scale_fill_manual(values=c('#91C4F2','#253C78'))
plot_grid(h0,featbox,ncol=2)
