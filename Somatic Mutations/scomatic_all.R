setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

#immune.combined<-subset(immune.combined,subset=Manuscript_Identity=='Dataset 1')

input.parent.dir <- "/home/rd796/scratch_public/Backup/Ruben/BaylorFASTQoutput/sample_out"
soup.batch.names <- list.files(file.path(input.parent.dir))
soup.batch.names <- soup.batch.names[!grepl('C',soup.batch.names)] #Remove IPF Cell Atlas
soup.batch.dir.paths <- file.path(input.parent.dir, soup.batch.names)

load('emptytbl.RData')
#load('emptytbl2.RData') #AT2S/AT2B
calltbl=fulltbl
    for(j in 1:length(soup.batch.names)){
    allmuts<-read.table(file=paste0(soup.batch.dir.paths[j],"/basecall/",soup.batch.names[j],".calling.step2.tsv"), sep='\t')
    allcall<-read.table(file=paste0(soup.batch.dir.paths[j],"/callable/",soup.batch.names[j],".coverage_cell_count.report.tsv"), header=T, sep='\t')
    allmuts<-allmuts[allmuts$V6=='PASS',]
    allmuts$V13<-as.numeric(allmuts$V13)
    allmuttbl<-allmuts %>% group_by(V7) %>% summarize(tot=sum(V13),num=length(V13))
    mtc<-match(allmuttbl$V7,allcall$Cell_types)
    allmuttbl<-cbind(allmuttbl,calls=allcall$DP[mtc+1])
    #allmuttbl<-cbind(allmuttbl,burd=allmuttbl$num/allmuttbl$calls)
    #allmuttbl<-allmuttbl[!(allmuttbl$V7 %in% c('T','B','DC','NK')),]
    allmuttbl$V7[allmuttbl$V7=='_Macrophage']='Alv. Macrophage'
    mtc2<-match(allmuttbl$V7,rownames(fulltbl))
    fulltbl[mtc2,j]<-allmuttbl$num
    calltbl[mtc2,j]<-allmuttbl$calls
    names(fulltbl)[j]<-soup.batch.names[j]
    names(calltbl)[j]<-soup.batch.names[j]
    fulltbl=cbind(fulltbl,rep(0,nrow(fulltbl)))
    calltbl=cbind(calltbl,rep(0,nrow(calltbl)))
    }
fulltbl<-fulltbl[,-(j+1)]
calltbl<-calltbl[,-(j+1)]
fulltbl[fulltbl==0]<-NA
calltbl[calltbl==0]<-NA
burdtbl=fulltbl/calltbl

#EXPORT GLOBAL DATA
load('imm_ages.RData')
ages=as.numeric(imm_ages$age[match(soup.batch.names,imm_ages$id)])
ages2=factor(cut(ages,c(0,40,60,100),c('young','mid','old')),levels=c('young','mid','old'))
gfg_data <- data.frame(ages=ages,mut=colSums(fulltbl,na.rm=T),call=colSums(calltbl,na.rm=T),burden=colSums(fulltbl,na.rm=T)/colSums(calltbl,na.rm=T))
#save(gfg_data,file='mutmatv2.RData')
#save(gfg_data,file='mutmat2v2.RData')
gfg_sample<-gfg_data

#sample versus age; not great, suffers from cell type bias
ggplot(gfg_data,aes(x=ages,y=log10(burden)))+geom_point()+theme_minimal()+labs(x='Age',y='log Mutation Burden')

#EXPORT CELL DATA
gfg_data <- data.frame(t(burdtbl))
#save(gfg_data,file='cellmutmatv2.RData')
#save(gfg_data,file='cellmutmat2v2.RData')

#GENERATE SUMMARY DATA
gfg_data=pivot_longer(gfg_data,colnames(gfg_data))
gfg_sum<-gfg_data %>% group_by(name) %>% summarise(avg=mean(value,na.rm=T),len=sum(!is.na(value)))
gfg_sum$name[gfg_sum$name=='Alv..Macrophage']='Alv. Macrophage'
gfg_sum=gfg_sum[match(rownames(burdtbl),gfg_sum$name),]
#OR...
rs<-rowSums(fulltbl,na.rm=T)/rowSums(calltbl,na.rm=T)
gfg_sum<-gfg_sum %>% mutate(totmut=rs)

#Plot per cell type
gfg_data <- data.frame((t(burdtbl)))
#gfg_data <- gfg_data[,-match(gfg_sum$name[gfg_sum$len<5],colnames(gfg_data))]
gfg_data <- gfg_data[,c('Aerocyte','Alv..Macrophage','AT1','AT2','Ciliated','gCap','Macrophage','Monocyte','NK','T')]
gfg_data=pivot_longer(gfg_data,colnames(gfg_data))
gfg_data=gfg_data[!is.na(gfg_data$value),]
gfg_data$name[gfg_data$name=='Alv..Macrophage']='Alv. Macrophage'
lvl=gfg_data %>% group_by(name) %>% summarise(mean=(mean((value))*1000000),med=median(value))
gfg_data$name=factor(gfg_data$name,levels=lvl$name[order(lvl$mean)])
ggplot(gfg_data,aes(x=name,y=(value)))+geom_boxplot(fatten=NULL)+geom_point()+coord_flip()+theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
labs(x='Cell Type',y='log Mutation Burden')+stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")

#NEW VERSION 5/7/25
gfg_data <- data.frame((t(burdtbl)))
#gfg_data <- gfg_data[,-match(gfg_sum$name[gfg_sum$len<5],colnames(gfg_data))]
gfg_data <- gfg_data[,c('Aerocyte','Alv..Macrophage','AT1','AT2','Ciliated','gCap','Macrophage','Monocyte','NK','T')]
gfg_data=pivot_longer(gfg_data,colnames(gfg_data))
gfg_data=gfg_data[!is.na(gfg_data$value),]
gfg_data$name[gfg_data$name=='Alv..Macrophage']='Alv. Macrophage'
lvl=gfg_data %>% group_by(name) %>% summarise(mean=(mean((value))*1000000),med=median(value)*1000000)
gfg_data$name=factor(gfg_data$name,levels=lvl$name[order(lvl$med)])
ggplot(gfg_data,aes(x=name,y=log10(value)))+geom_boxplot()+geom_point()+coord_flip()+ylim(-6.25,-3.75)+
theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(x='Cell Type',y='log Mutation Burden')

#per sample per cell type
gfg_data <- data.frame(Age=ages,log10(t(burdtbl)))
gfg_data=gfg_data[,c(T,!gfg_sum$len<3)] #deleting row #1
plist=lapply(2:11,function(x){ggplot(gfg_data,aes(x=Age,y=eval(as.name(paste(colnames(gfg_data)[x])))))+geom_point()+theme_minimal()})
library(gridExtra)
do.call("grid.arrange", c(plist, ncol=5))

#sample-cell type basis
load('fill_df.RData')
fill_df$predicted.id[fill_df$predicted.id=='Alv. Fibroblast']='Fibroblast'
gfg_data <- data.frame(Age=ages,log10(t(burdtbl)))
#gfg_data <- gfg_data[,-match(gfg_sum$name[gfg_sum$len<5],colnames(gfg_data))]
gfg_data <- gfg_data[,c('Age','Aerocyte','Alv..Macrophage','AT1','AT2','Ciliated','gCap','Macrophage','Monocyte')]
gfg_data <- gfg_data[gfg_sample$mut<300,] #READ: REMOVE 367 d/t high burden
#gfg_data <- gfg_data[gfg_sample$mut<2000,] #READ: REMOVE 367 d/t high burden
apply(gfg_data,2,function(x){cor(x,gfg_data$Age,use='pairwise.complete.obs')})
#apply(gfg_data,2,function(x){corPvalueStudent(cor(x,gfg_data$Age,use='pairwise.complete.obs'),sum(!is.na(x)))})
gfg_data <- pivot_longer(gfg_data,!Age)
gfg_data <- gfg_data[!is.na(gfg_data$value),]
gfg_data$name[gfg_data$name=='Alv..Macrophage']='Alv. Macrophage'
gfg_data$name=factor(gfg_data$name,levels=sort(unique(gfg_data$name)))
gfg_data <- gfg_data %>% mutate(color1=fill_df$color[match(gfg_data$name,fill_df$predicted.id)])
mycolors=unlist(fill_df$color)
names(mycolors)=fill_df$predicted.id
colfun <- scale_colour_manual(name='Cell Type',values=mycolors)
#p2<-ggplot(gfg_data,aes(x=Age,y=value,color=name))+geom_point()+theme_minimal()+scale_colour_manual(values=fill_df$color[match(gfg_data$name,fill_df$predicted.id)])+ylab('Mutation Burden')+labs(color='Cell Type')
#p2<-ggplot(gfg_data,aes(x=Age,y=value,color=name))+geom_point()+theme_minimal()+scale_colour_manual(values=fill_df$color[match(gfg_data$name,fill_df$predicted.id)])+ylab('Mutation Burden')+labs(color='Cell Type')+geom_smooth(method='lm',se=F) 
p2<-ggplot(gfg_data,aes(x=Age,y=value,color=name))+geom_point()+geom_smooth(method='gam',formula = y ~ s(x,k=4),se=F) +theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+colfun+ylab('log Mutation Burden')+labs(color='Cell Type')
p2
cor(gfg_data$Age,gfg_data$value,use='pairwise.complete.obs')
library(glmmTMB)
  f2 <- glmmTMB(value ~ Age +  (1 | name), data = gfg_data)
summary(f2)

###########
# SENMAYO #
###########
tbl=data.frame(sen=immune.combined$SenMayo1,cell=immune.combined$predicted.id,id=immune.combined$orig.ident)
sen_ct=tbl %>% group_by(cell) %>% summarise(sen_ct=mean(sen))
sen_id=tbl %>% group_by(id) %>% summarise(sen_id=mean(sen))

#per sample
gfg_data <- data.frame(ages=ages,call=colSums(calltbl,na.rm=T),mut=colSums(fulltbl,na.rm=T),sen=sen_id$sen_id)
gfg_data=cbind(gfg_data,b=gfg_data$mut/gfg_data$call)
gfg_data=gfg_data[!gfg_sum$len<2,]
ggplot(gfg_data,aes(x=sen,y=b))+geom_point()+theme_minimal()

#per cell type
library(ggrepel)
gfg_data <- data.frame(call=rowSums(calltbl,na.rm=T),mut=rowSums(fulltbl,na.rm=T))
gfg_data=gfg_data %>% mutate(b=mut/call,sen=sen_ct$sen_ct[match(rownames(gfg_data),sen_ct$cell)])
gfg_data=gfg_data[!gfg_sum$len<2,]
ggplot(gfg_data,aes(x=sen,y=b))+geom_point()+geom_label_repel(aes(label=rownames(gfg_data)),size=3)+theme_minimal()+scale_y_continuous(limits = c(0, 0.002))+
labs(x='Senescence Score',y='Mutation Burden')
