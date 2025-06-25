setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

input.parent.dir <- "/home/rd796/scratch_public/Backup/Ruben/BaylorFASTQoutput/sample_out"
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
colnames(matlistsum2)[2]<-

ggplot(matlistsum2,aes(x=source,y=Count,fill=Func.refGene))+geom_bar(position="fill", stat="identity")+facet_wrap(~Age,nrow=1,scales='free_x')+
theme_minimal()+theme(axis.text.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(x='Sample',y='Proportion')
