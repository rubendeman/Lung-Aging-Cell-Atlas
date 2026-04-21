library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(gprofiler2)

#Add path to GLMM output
path=?

#read glmmTMB output
comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-as.list(unlist(comp)); names(cell.types)<-cell.types
datal=lapply(cell.types,function(x){read.table(paste0(path,x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(cell.types,function(x){read.table(paste0(path,x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(datal,FUN = setNames,nm=c('Beta0','Beta1','p_val','p_val_adj'))
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})

names(datal)<-str_replace_all(names(datal), c("Alv_Macrophage" = "Alv. Macrophage", "Adv_Fibroblast" = "Adventitial Fibroblast", "Alv_Fibroblast"="Alv. Fibroblast"))
datalsig<-lapply(datal,function(x){x[x$p_val<0.05,]})
datalfdr<-lapply(datal,function(x){x[x$p_val_adj<0.05,]})

#FIGURE 2
plotdata<-data.frame(ID=rep(names(datal),sapply(datal,nrow)),do.call('rbind',datal))
plotdata=plotdata %>% mutate(Significant=p_val_adj<0.05)
plotdata=plotdata[plotdata$p_val<0.05,]

load('fill_df.RData')
ups_data=lapply(datal,function(x){rownames(x[x$p_val<0.05,])})
lmat=data.frame(Cell=names(ups_data),Diff=sapply(ups_data,length),Number=sapply(datal,nrow))
lmat=lmat %>% mutate(PerDiff=Diff/Number)

plotdata$Beta1=as.numeric(plotdata$Beta1)*-1 #Make positive
g1<-ggplot(plotdata %>% arrange(Significant), aes(x=factor(ID,levels=lmat$Cell[order(lmat$PerDiff)]), y=Beta1,colour=Significant))+geom_point(position=position_jitter())+
scale_colour_manual(values=c('grey','blue'))+scale_y_continuous(limits=c(-5,5))+scale_x_discrete(position='top')+
theme_minimal()+theme(axis.line.x=element_line(size=1),axis.text.y.right = element_text(color='black',hjust=0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+coord_flip()+NoLegend()+labs(x=NULL)
lmat=lmat %>% mutate(color=fill_df$color[match(lmat$Cell,fill_df$predicted.id)])
lmat=lmat[order(lmat$PerDiff,decreasing=T),]
lmat$Cell=factor(lmat$Cell,levels=lmat$Cell)
lmat$PerDiff=lmat$PerDiff*100
g2<-ggplot(lmat, aes(y = PerDiff, x = factor(Cell, levels=Cell[order(PerDiff)]), fill=Cell)) +
    geom_bar(stat = "identity")+theme_minimal()+theme(axis.line.x=element_line(size=1),axis.text.y =element_blank(), axis.title.y = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+NoLegend()+
    labs(x="Cell Type", y = "% Differentially Expressed Genes", fill=NULL)+coord_flip()+scale_fill_manual(values=lmat$color)
g1+g2

#UPSETS
library(UpSetR)
data_upset=lapply(datal,function(x){x$Gene[x$p_val<0.05]}) #UP

#Global UPSET
upset(fromList(data_upset),nsets=length(data_upset),order.by='freq')
tmp=unlist(data_upset)
tmp=as.data.frame(table(tmp))
feats=as.character(tmp$tmp[tmp$Freq>4])
noquote(feats)

#Endothelial UPSET
ct=comp[['Endothelial']]
u1<-upset(fromList(data_upset[ct]),order.by='freq',empty.intersections="on",text.scale=1.2)
u1

#Capillary Overlap
ct=c('gCap','Aerocyte')

data_upset=lapply(datal,function(x){x$Gene[x$p_val<0.05&x$Beta1<0]}) #UP
upgenes=Reduce(intersect,data_upset[ct])
rnk=order(colMeans(rbind(match(upgenes, data_upset[[ct[1]]]),match(upgenes, data_upset[[ct[2]]]))))

data_upset=lapply(datal,function(x){x$Gene[x$p_val<0.05&x$Beta1>0]}) #DN
dngenes=Reduce(intersect,data_upset[ct])
rnk=order(colMeans(rbind(match(dngenes, data_upset[[ct[1]]]),match(dngenes, data_upset[[ct[2]]]))))

#GO
golist=upgenes
go_out<-gost(query=golist, organism = "hsapiens")
kres=go_out$result[go_out$result$source%in%c('KEGG','REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>10){kmat<-kmat[1:12,]}
kmat=kmat[1:5,]
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
u1<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) + geom_col() +
    scale_fill_gradient(low='grey',high='yellow',limits=c(0,max(-log10(kmat$p_value))))+ geom_text(aes(label = term_name,hjust=1)) + 
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='Enrichment Term')+ guides(fill=guide_legend(title="-log10(p-value)"))
u1

golist=dngenes
go_out<-gost(query=golist, organism = "hsapiens")
kres=go_out$result[go_out$result$source%in%c('KEGG','REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>10){kmat<-kmat[1:12,]}
kmat=kmat[c(1,2,8,10,11),]
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
kmat=kmat %>% mutate(place = if_else(row_number() <4, 1, 0))
kmat=kmat %>% mutate(color = if_else(row_number() <4, 'white', 'black')) 
u2<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) + geom_col() +
    scale_fill_gradient(low='grey',high='purple',limits=c(0,max(-log10(kmat$p_value))))+ geom_text(aes(label = term_name),hjust=kmat$place,colour=kmat$color) + 
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='Enrichment Term')+ guides(fill=guide_legend(title="-log10(p-value)"))
u2
plot_grid(u1,u2,ncol=1)
