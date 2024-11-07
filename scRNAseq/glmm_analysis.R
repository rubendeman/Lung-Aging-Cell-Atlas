setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(gprofiler2)

#read glmmTMB output
comp<-list(Epithelial=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),Endothelial=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),Mesenchymal=c('Adv_Fibroblast','Alv_Fibroblast','Myofibroblast','SMC','Pericyte'),Myeloid=c('Monocyte','Macrophage','Alv_Macrophage'),Lymphoid=c('B','T','Mast','DC','NK'))
cell.types<-as.list(unlist(comp)); names(cell.types)<-cell.types
#datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/project/ageproj/HPC_GLMM_AGE_GENES/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
#datal=lapply(datal,FUN = setNames,nm=c('Beta0','Beta1','p_val','p_val_adj'))
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})

names(datal)<-str_replace_all(names(datal), c("Alv_Macrophage" = "Alv. Macrophage", "Adv_Fibroblast" = "Adventitial Fibroblast", "Alv_Fibroblast"="Alv. Fibroblast"))

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
#pdf('fig2.pdf',width=18,height=12) #2X
#plot_grid(g1,g2,g3,g4,ncol=2)
#dev.off()

#DGE
bla<-datal[['gCap']]
bladn<-bla[bla$Beta1<0,]
blaup<-bla[bla$Beta1>0,]
upgenes<-bladn$Gene[1:25]
dngenes<-blaup$Gene[1:25]
noquote(bladn$Gene[1:100])
noquote(blaup$Gene[1:100])
#sapply(datal,function(x){x$p_val_adj[x$Gene=='FBXW7']})

golist <- bladn$Gene[bladn$p_val_adj<0.05]
golist <- blaup$Gene[blaup$p_val_adj<0.05]

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
u1<-upset(fromList(data_upset[ct]),order.by='freq',empty.intersections="on")
u1

#Capillary Overlap
ct=c('gCap','Aerocyte')

data_upset=lapply(datal,function(x){x$Gene[x$p_val<0.05&x$Beta1<0]}) #UP
upgenes=Reduce(intersect,data_upset[ct])
rnk=order(colMeans(rbind(match(upgenes, data_upset[[ct[1]]]),match(upgenes, data_upset[[ct[2]]]))))
upgenes=upgenes[rnk][1:21]
data_upset=lapply(datal,function(x){x$Gene[x$p_val<0.05&x$Beta1>0]}) #DN
dngenes=Reduce(intersect,data_upset[ct])
rnk=order(colMeans(rbind(match(dngenes, data_upset[[ct[1]]]),match(dngenes, data_upset[[ct[2]]]))))
dngenes=dngenes[rnk][1:21]

#Query autophagy genes
ct=c('Aerocyte','gCap','Arterial','Venous','AT1','AT2')
library(openxlsx)
ubq=lapply(list(2,3,4,5,6),function(x){read.xlsx('Ge_UBQ.xlsx',sheet=x,cols=1,skipEmptyRows=T)})
names(ubq)<-c('E1','E2','E3 Activity','E3 Adapter','DUB')
auto=lapply(list(1,2,3,4,5,6,7),function(x){read.xlsx('Bordi_Autophagy.xlsx',sheet=x,cols=1,skipEmptyRows=T)})
names(auto)<-c('mTOR','Core','Reg','Mitophagy','Docking','Lysosome','Lysosome Genes')
auto=c(ubq,auto)
auto2=lapply(ct,function(y){lapply(auto,function(x){na.omit(data.frame(gene=datal[[y]]$Gene[match(x[,1],datal[[y]]$Gene)],beta1=datal[[y]]$Beta1[match(x[,1],datal[[y]]$Gene)]*-1,pv=datal[[y]]$p_val[match(x[,1],datal[[y]]$Gene)]))})})
names(auto2)<-ct
auto2[["gCap"]][["mTOR"]]

#Query marker genes
ct=c('Aerocyte','gCap','Arterial','Venous','AT1','AT2')
auto=c(Cap=list(data.frame(Gene=c('IL7R','FCN3','EDN1','APLNR','TEK','CA4','EDNRB','FENDRR','VIPR1','APLNR','APLN','HPGD'))),Art=list(data.frame(Gene=c('DKK2','GJA5','ACKR1'))),Pan=list(data.frame(Gene=c('PECAM1','SOX17','CD34','CDH5','VWF'))),
AT1=list(data.frame(Gene=c('HOPX','AGER','RTKN2'))),AT2=list(data.frame(Gene=c('SFTPC','LAMP3','SLC34A2','SFTPB','SFTPA1'))))
auto2=lapply(ct,function(y){lapply(auto,function(x){na.omit(data.frame(gene=datal[[y]]$Gene[match(x[,1],datal[[y]]$Gene)],beta1=datal[[y]]$Beta1[match(x[,1],datal[[y]]$Gene)],pv=datal[[y]]$p_val[match(x[,1],datal[[y]]$Gene)]))})})
names(auto2)<-ct
auto2[["AT1"]][["AT1"]]

#Query NK collab
ct=c('AT1','AT2','Alv. Fibroblast','Adventitial Fibroblast','SMC','Myofibroblast')
auto=c(genes=list(data.frame(Gene=c('MIA3','MIA2','SEC23A','SEC23B','SEC13','SEC31A','SAR1A','SAR1B','SEC23IP','TFG','TUG1','ALG2','SERPINH1','SEC16A','SEC16B','PREB',
'RINT1','ZW10','NBAS','SCFD1','VPS13B','ATG2A','SEC22A','SEC22B','GOLGA6A','DYM','TXNDC5','EXPH5','CUL3','KLHL12','TMEM131','TMEM39A','TMED1','TMED2','TMED3','TMED4','TMED5','TMED6','TMED7','TMED9','TMED10','SURF4',
'CNI1','TRAPPC3','TRAPPC2','RAB1A','RAB1B','RAB2A','USO1'))))
auto2=lapply(ct,function(y){lapply(auto,function(x){(data.frame(gene=datal[[y]]$Gene[match(x[,1],datal[[y]]$Gene)],beta1=datal[[y]]$Beta1[match(x[,1],datal[[y]]$Gene)],pv=datal[[y]]$p_val[match(x[,1],datal[[y]]$Gene)]))})})
names(auto2)<-ct
auto3=as.data.frame(auto2)

#GO
golist=upgenes
go_out<-gost(query=golist, organism = "hsapiens");
kres=go_out$result[go_out$result$source%in%c('KEGG','REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>10){kmat<-kmat[1:12,]}
kmat=kmat[1:6,]
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
u2<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) + geom_col() +
    scale_fill_gradient(low='grey',high='yellow',limits=c(0,max(-log10(kmat$p_value))))+ geom_text(aes(label = term_name,hjust=1)) + 
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='REAC Term')+ guides(fill=guide_legend(title="-log10(p-value)"))
u2

golist=dngenes
go_out<-gost(query=golist, organism = "hsapiens");
kres=go_out$result[go_out$result$source%in%c('KEGG','REAC'),]
kres=kres[order(kres$p_value,decreasing=F),]
kmat=kres[,c(3,11)]
if(nrow(kmat)>10){kmat<-kmat[1:12,]}
kmat$term_name[2]='Respiratory electron transport, ATP synthesis'
kmat$term_name[3]='Chemical carcinogenesis - ROS'
kmat$term_name[7]='TCA cycle and respiratory electron transport'
kmat=kmat[c(1,2,3,6,7,12),]
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
u3<-ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) + geom_col() +
    scale_fill_gradient(low='grey',high='purple',limits=c(0,max(-log10(kmat$p_value))))+ geom_text(aes(label = term_name,hjust=1)) + 
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='REAC Term')+ guides(fill=guide_legend(title="-log10(p-value)"))
u3
plot_grid(u2,u3,ncol=1)
