setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(stats)
library(stringr)
library(WGCNA)
library(ggpubr)
load('evennewerstep1.RData') #load bulk

################################
# Boxplot - alongside scRNAseq #
################################
load('prop12_15.RData')
comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv. Macrophage','DC'),Lymphoid=c('B','T','Mast','NK'))
cell.types<-unlist(comp)

ind=tempprop$Est.prop.weighted
demo_sub<-as.data.frame(lapply(comp,function(x){rowSums(ind[,x],na.rm=T)}))
clinprop=clindata5

#for (i in 1:5){ind[,comp[[i]]]<-ind[,comp[[i]]]/demo_sub[,i]} #Convert ind from raw counts to sub props
#id_cond<-ind[,2]>0.05&ind[,16]>0
#idsub=rownames(ind)[id_cond]
#ind=ind[id_cond,]
#clinprop=clinprop[id_cond,]

p=list()
for (i in 1:5){
cormat2<-na.omit(as.data.frame(pivot_longer(as.data.frame(cbind(age=clinprop$AGE,ind[,comp[[i]]])),!age)))
vage<-cormat2$name
vage[cormat2$age<61]<-'Young'
vage[cormat2$age>60]<-'Aged'
cormat2$age<-vage
if(i==1){epidata=cormat2}
p[[i]]<-ggplot(cormat2, aes(y = value, x = name,fill=factor(age,levels=c('Young','Aged')))) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(size=12,angle = 45, vjust = 1, hjust=1)) +
    labs(x=NULL, y = "Gene Expression")+guides(fill=guide_legend(title=NULL))+geom_point(position=position_jitterdodge(jitter.width=0.1))+scale_fill_manual(values=c('#91C4F2','#253C78'))#+stat_compare_means(aes(label=after_stat(p.signif)))
}
plot_grid(plotlist=p) +theme(plot.margin = margin(1,1,1,1, "cm"))

wilcox.test(epidata$value[epidata$name=='AT2'&epidata$age=='Aged'],epidata$value[epidata$name=='AT2'&epidata$age=='Young'],alternative='less')
