setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

#Load Seurat
load('agingseurat.RData')
#immune.combined$predicted.id[immune.combined$predicted.id=='AT2B']<-'AT2'
#immune.combined$predicted.id[immune.combined$predicted.id=='AT2S']<-'AT2'
load('fill_df.RData')

comp<-list(Epithelial=c('AT1','AT2B','AT2S','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage','DC'),Lymphoid=c('B','T','Mast','NK'))
cell.types<-unlist(comp)

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
id_cond<-rowSums(cbind(ind[,2],ind[,3]),na.rm=T)>5
idsub=rownames(ind)[id_cond]
ind=ind[id_cond,]
demo_sub=demo_sub[id_cond,]

demo=as.data.frame(cbind(count=apply(ind,1,function(x){sum(x,na.rm=T)}),age=paste(immune.combined@meta.data$agebin,immune.combined@meta.data$smoke)[match(idsub,immune.combined@meta.data$orig.ident)]))

#ind<-ind/demo$count #Convert ind from raw counts to props
#or
for (i in 1:5){ind[,comp[[i]]]<-ind[,comp[[i]]]/demo_sub[,i]} #Convert ind from raw counts to sub props

p=list()
for (i in 1:5){
cormat2<-na.omit(as.data.frame(pivot_longer(as.data.frame(cbind(age=demo$age,ind[,comp[[i]]])),!age)))
cormat2$value=as.numeric(cormat2$value)
if(i==1){epidata=cormat2}
p[[i]]<-ggplot(cormat2, aes(y = value, x = name,fill=age)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Cell Type Proportion")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))+stat_compare_means(method='wilcox.test')
}
plot_grid(plotlist=p)+theme(plot.margin = margin(1,1,1,1, "cm"))

wilcox.test(epidata$value[epidata$name=='AT2'&epidata$age=='Aged'],epidata$value[epidata$name=='AT2'&epidata$age=='Young'],alternative='less')

wilcox.test(cormat2$value[cormat2$name=='SMC'&cormat2$age=='Yes'],cormat2$value[cormat2$name=='SMC'&cormat2$age=='No'],alternative='less')
