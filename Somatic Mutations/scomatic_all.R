setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)

immune.combined<-subset(immune.combined,subset=Manuscript_Identity=='Dataset 1')

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
    allmuts$V13<-as.numeric(allmuts$V13) #OR 11/13
    allmuttbl<-allmuts %>% group_by(V7) %>% summarize(tot=sum(V13)) #OR 11/13
    mtc<-match(allmuttbl$V7,allcall$Cell_types)
    allmuttbl<-cbind(allmuttbl,calls=allcall$DP[mtc+1])

    #allmuttbl$V7[allmuttbl$V7%in%c('Fibroblast','Myofibroblast','SMC')]='Fibroblast' #***************
    #allmuttbl$V7[allmuttbl$V7%in%c('Arterial','Venous','Peribronchial')]='Arterial' #***************
    #allmuttbl = allmuttbl %>% group_by(V7) %>% summarise_each(list(sum))

    allmuttbl<-cbind(allmuttbl,burd=allmuttbl$tot/allmuttbl$calls)
    allmuttbl<-allmuttbl[!(allmuttbl$V7 %in% c('T','B','DC','NK')),]
    allmuttbl$V7[allmuttbl$V7=='_Macrophage']='Alv. Macrophage'
    mtc2<-match(allmuttbl$V7,rownames(fulltbl))
    fulltbl[mtc2,j]<-allmuttbl$tot
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

load('imm_ages.RData')
ages=as.numeric(imm_ages$age[match(soup.batch.names,imm_ages$id)])
ages2=factor(cut(ages,c(0,40,60,100),c('young','mid','old')),levels=c('young','mid','old'))

#sample versus age; not great, suffers from cell type bias
gfg_data <- data.frame(ages=ages,call=colSums(calltbl,na.rm=T),mut=colSums(fulltbl,na.rm=T),b=colSums(fulltbl,na.rm=T)/colSums(calltbl,na.rm=T))
ggplot(gfg_data,aes(x=ages,y=log10(b)))+geom_point()+theme_minimal()+labs(x='Age',y='log Mutation Burden')

#SUMMARY TABLE
gfg_data <- data.frame(t(burdtbl))
gfg_data=pivot_longer(gfg_data,colnames(gfg_data))
gfg_sum<-gfg_data %>% group_by(name) %>% summarise(avg=mean(value,na.rm=T),len=sum(!is.na(value)))
gfg_sum$name[gfg_sum$name=='Alv..Macrophage']='Alv. Macrophage'
gfg_sum=gfg_sum[match(rownames(burdtbl),gfg_sum$name),]

#per cell type
gfg_data <- data.frame(log10(t(burdtbl)))
#gfg_data <- gfg_data[,-match(gfg_sum$name[gfg_sum$len<5],colnames(gfg_data))]
gfg_data <- gfg_data[,c('Aerocyte','Alv..Macrophage','AT1','AT2','Ciliated','gCap','Macrophage','Monocyte')]
gfg_data=pivot_longer(gfg_data,colnames(gfg_data))
gfg_data=gfg_data[!is.na(gfg_data$value),]
gfg_data$name[gfg_data$name=='Alv..Macrophage']='Alv. Macrophage'
lvl=gfg_data %>% group_by(name) %>% summarise(mean=(mean((value))*1000000))
gfg_data$name=factor(gfg_data$name,levels=lvl$name[order(lvl$mean)])
ggplot(gfg_data,aes(x=name,y=(value)))+geom_boxplot(fatten=NULL)+geom_point()+coord_flip()+theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
labs(x='Cell Type',y='log Mutation Burden')+stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")

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
gfg_data <- gfg_data[,-match(gfg_sum$name[gfg_sum$len<5],colnames(gfg_data))]
gfg_data <- gfg_data[,c('Age','Aerocyte','Alv..Macrophage','AT1','AT2','Ciliated','gCap','Macrophage','Monocyte')]
gfg_data <- gfg_data[gfg_sample$mut<2000,]
apply(gfg_data,2,function(x){cor(x,gfg_data$Age,use='pairwise.complete.obs')})
apply(gfg_data,2,function(x){corPvalueStudent(cor(x,gfg_data$Age,use='pairwise.complete.obs'),sum(!is.na(x)))})
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

################
# READ GLMMTMB #
################
#read glmmTMB output
datal=read.table(paste0('/home/rd796/palmer_scratch/global/age-glmmTMB_global_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)
datal=datal[!(is.na(datal$Beta1)|is.na(datal$p_val)|!is.na(datal$blank)),]
datal=mutate(datal,p_val_adj=p.adjust(datal$p_val,method='fdr',n=nrow(datal)))
datal=datal[order(datal$p_val,decreasing=F),]
noquote(datal$Gene[datal$Beta1<0][1:40])
noquote(datal$Gene[datal$Beta1<0&datal$p_val<0.05])

##########
# GLOBAL #
##########
gfg_data <- data.frame(ages=ages,call=colSums(calltbl,na.rm=T),mut=colSums(fulltbl,na.rm=T))
gfg_data=cbind(gfg_data,b=gfg_data$mut/gfg_data$call)
ann <- gfg_data$b
rnames=rownames(gfg_data)
tmp=!is.na(ann)
ann<-ann[tmp]
ann=factor(cut(ann,c(0,median(ann),100),c('Low Mutation','High Mutation')),levels=c('Low Mutation','High Mutation')) #OPTIONAL; BINARY MUT
rnames<-rnames[tmp]
ann=data.frame('Mutation Burden'=(ann)) #log10
rownames(ann)=rnames
ann=ann %>% mutate(Age=immune.combined$age[match(rownames(ann),immune.combined$orig.ident)])

small=subset(immune.combined,subset=orig.ident%in%rnames)
small=AddMetaData(small,metadata=ann$'Mutation.Burden'[match(small$orig.ident,rnames)],col.name='mut')

VlnPlot(small,features='RAD50',split.by='orig.ident',group.by='mut')

#CDKN2A Plot
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features='CDKN2A'))
datalist=data.frame(exp=t(avg),mut=ann$Mutation.Burden)
ggplot(datalist, aes(x=mut, y=exp))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+theme_classic()+theme(text = element_text(size=15))+
labs(x=NULL, y = "CDKN2A Expression")+scale_y_continuous(limits=c(0,0.35))+stat_compare_means(label='p.signif')

feats=c(datal$Gene[datal$Beta1<0][1:20],datal$Gene[datal$Beta1>0][1:20])
rowsplit=data.frame(Gene=rep(c('Decreased','Increased'),each=length(feats)/2))
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features=feats))

################
# MAKE HEATMAP #
################
ann=ann[(rownames(ann) %in% substring(colnames(avg),5)),,F] #sanity
names(ann)[1]<-'Mutations'
hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=circlize::colorRamp2(c(0, 100), c("white", "blue")))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann, col=hcols,which='col',gap = unit(1, 'mm'))

rowfeats=c('COX4I1','NDUFB8','ATP','RAD50','PRKDC','ANAPC1','SEL1L','USO1','RAB3GAP1','COQ10')
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avg) %in% rowfeats), labels = rownames(avg)[rownames(avg)%in%rowfeats]))
ha=NULL
avg2=t(apply(avg,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht=Heatmap(na.omit(avg2), name = "Expression",  column_split=data.frame(factor(ann$'Mutations')),cluster_columns = F, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F,cluster_row_slices = F, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = T,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=T,column_title = NULL, row_split=rowsplit,use_raster = F)
draw(ht, padding = unit(c(2, 2, 8, 4), "mm"),merge_legend=T)

################
# BY CELL TYPE #
################
#read glmmTMB output
cell.types<-c('AT2','gCap')
datal=lapply(cell.types,function(x){read.table(paste0('/home/rd796/palmer_scratch/cell/age-glmmTMB_',x,'_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)})
datal=lapply(datal,function(x){x[!(is.na(x$Beta1)|is.na(x$p_val)|!is.na(x$blank)),]})
datal=lapply(datal,function(x){mutate(x,p_val_adj=p.adjust(x$p_val,method='fdr',n=nrow(x)))})
datal=lapply(datal,function(x){x[order(x$p_val,decreasing=F),]})
names(datal)<-cell.types

gfg_data <- data.frame(t(burdtbl))

ht=NULL
for (i in cell.types){
ann <- gfg_data[,i]
rnames=rownames(gfg_data)
tmp=!is.na(ann)
if(i=='gCap'){tmp[18]=F}
ann<-ann[tmp]
ann=factor(cut(ann,c(0,median(ann),100),c('Low Mutation','High Mutation')),levels=c('Low Mutation','High Mutation'))
rnames<-rnames[tmp]
ann=data.frame('Mutation Burden'=(ann)) #log10
rownames(ann)=rnames
if(i=='gCap'){ann$Mutation.Burden[rownames(ann)=="460"]="High Mutation"}
ann=ann %>% mutate(Age=immune.combined$age[match(rownames(ann),immune.combined$orig.ident)],'Cell Type'=i)

small<-subset(immune.combined,subset=predicted.id==i&orig.ident %in% rnames)
small=AddMetaData(small,metadata=ann$'Mutation.Burden'[match(small$orig.ident,rnames)],col.name='mut')

VlnPlot(small,features='MT-CO1',split.by='orig.ident',group.by='mut')

feats=c(datal[[i]]$Gene[datal[[i]]$Beta1>0][1:20],datal[[i]]$Gene[datal[[i]]$Beta1<0][1:20])
noquote(feats)
rowfeats=c('VIPR1','MACROD1','FENDRR')
rowsplit=data.frame(Gene=rep(c('Increased','Decreased'),each=length(feats)/2))
avg=as.data.frame(AverageExpression(small,group.by='orig.ident',assays='RNA',features=feats))

################
# MAKE HEATMAP #
################
ann=ann[(rownames(ann) %in% substring(colnames(avg),5)),,F] #sanity
names(ann)[1]<-'Mutations'
hcols=list('Cell Type'=fill_df$color[match(i,fill_df$predicted.id)],Age=circlize::colorRamp2(c(0, 100), c("white", "blue")),
'Mutations'=c('High Mutation'='#A20021','Low Mutation'='#F52F57'))
names(hcols$'Cell Type')=i
colAnn <- HeatmapAnnotation(df = ann, col=hcols,which='col',gap = unit(1, 'mm'))
ha = rowAnnotation(foo = anno_mark(at = which(rownames(avg) %in% rowfeats), labels = rownames(avg)[rownames(avg)%in%rowfeats]))
ha=NULL
avg2=t(apply(avg,1,function(x){(x-min(x))/(max(x)-min(x))}))
col_fun = circlize::colorRamp2(c(0, 1), c( "black", "#FFFF00"))
ht0=Heatmap(na.omit(avg2), name = "Expression",  column_split=data.frame(factor(ann$'Mutations')),cluster_columns = T, show_column_dend = FALSE,top_annotation=colAnn,right_annotation=ha,
            cluster_column_slices = F, column_title_gp = gpar(fontsize = 12), column_gap = unit(0.5, "mm"), cluster_rows = TRUE,
            show_row_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 12), row_title_rot = 0, column_title_rot = 45,
            show_column_names = FALSE, show_row_names=T,column_title = NULL, row_split=rowsplit,use_raster = F)
            ht=ht+ht0
}
draw(ht, padding = unit(c(2, 2, 8, 6), "mm"),auto_adjust=F,merge_legend=T)
