library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

input.parent.dir <- ?
soup.batch.names <- list.files(file.path(input.parent.dir))
soup.batch.names <- soup.batch.names[!grepl('C',soup.batch.names)] #Remove IPF Cell Atlas
soup.batch.dir.paths <- as.list(file.path(input.parent.dir, soup.batch.names))

matlist<-lapply(soup.batch.dir.paths,function(x){tryCatch(read.csv(file=paste0(x,"/annovar/sample.variants.annovar.hg38_multianno.csv")),error=function(e) NULL)})
names(matlist)<-soup.batch.names

matall<-do.call("rbind",matlist)
ggplot(matall,aes(x=Func.refGene))+geom_bar()+theme_minimal()+labs(x='Type',y='Count')

#With error bars for samples
matlistsum<-bind_rows(matlist[-23],.id="source")
matlistsum2<-matlistsum%>%group_by(source,Func.refGene)%>%summarise(Count=length(Func.refGene))

matlistsum2$Age=immune.combined$agebin[match(matlistsum2$source,immune.combined$orig.ident)]
#colnames(matlistsum2)[2]<-

ggplot(matlistsum2,aes(x=source,y=Count,fill=Func.refGene))+geom_bar(position="fill", stat="identity")+facet_wrap(~Age,nrow=1,scales='free_x')+
theme_minimal()+theme(axis.text.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(x='Sample',y='Proportion')

#Exonic types of mutations
matlistsum<-bind_rows(matlist[-23],.id="source")
matlistsum$Age=immune.combined$agebin[match(matlistsum$source,immune.combined$orig.ident)]
matlistsum<-matlistsum%>%group_by(Age,ExonicFunc.refGene)%>%summarise(Count=length(Func.refGene))
matlistsum=matlistsum[!matlistsum$ExonicFunc.refGene %in% c('.'),]

#Signatures
#golist=matall$Gene.refGene[matall$Func.refGene=='exonic']
golist=matall$Gene.refGene[matall$ExonicFunc.refGene=='nonsynonymous SNV']
library('gprofiler2')
go_out<-gost(query=golist, organism = "hsapiens")
kres=go_out$result[go_out$result$source=='KEGG',]
#kres=go_out$result[go_out$result$source=='KEGG'&go_out$result$term_size<200,]
kmat=kres[,c(3,11)]
if(nrow(kmat)>5){kmat<-kmat[1:5,]}
kmat$term_name<-factor(kmat$term_name,levels=rev(kmat$term_name))
kmat<-kmat[-4,]
kmat=kmat %>% mutate(place = if_else(row_number() <4, 1, 0))
kmat=kmat %>% mutate(color = if_else(row_number() <4, 'white', 'black')) 
ggplot(kmat, aes(x = -log10(p_value), y = term_name, fill=-log10(p_value))) +
    geom_col()+scale_fill_gradient(low="lightblue",high="blue",limits=c(0,max(-log10(kmat$p_value))))+
    geom_text(aes(label = term_name),hjust = kmat$place,colour=kmat$color) +
    theme_classic()+theme(legend.position="right",axis.text.y=element_blank())+labs(x='-log10(p-value)',y='KEGG Term')+ guides(fill=guide_legend(title="-log10(p-value)"))

noquote(golist)
