setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(ggplot2)
library(dplyr)
library(WGCNA)
library(DESeq2)
library(EnhancedVolcano)
load('evennewerstep1.RData') #load demographic data
load('loggedmes.RData') #load module eigengenes
load('taball.RData') #load bulk mutation data
load('prop12_15.RData') #load music props
head(tab_all,3)
tbl=tab_all[tab_all$tissue == "Lung",]
mutburden=as.data.frame(table(tbl$gtexIds))

comp<-list(epi=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),endo=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),mes=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),imm=c('Monocyte','Macrophage','Alv. Macrophage'),lymph=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)

#Load bulk
datExpr=read.csv('raw.csv')
ids=read.csv('idgenename.csv')
dp=duplicated(ids$Description)
datExpr=datExpr[!dp,]
datExpr=datExpr[,c('id',rownames(MEs))]
ids=ids[!dp,]
probes2annot = match(datExpr$id, ids$id)
rownames(datExpr) = ids$Description[probes2annot]
datExpr=datExpr[,-1]

#Plot
tbl[,6]=paste(tbl[,3], tbl[,4], sep=">")
#tbl[,6]=vapply(paste(tbl[,3], tbl[,4], sep=">"), function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=''), '')
tbl2=as.data.frame(table(tbl[,6]))
tbl2[c(1:3,7:9),1]=c('T>G','T>C','T>A','C>T','C>G','C>A')
tbl2=tbl2 %>% group_by(Var1) %>% summarise(Freq=sum(Freq))
ggplot(tbl2, aes(y = Freq, x = Var1, fill=Var1)) + geom_col() + labs(fill='Base Substitution')+theme_classic()+theme(line = element_blank(),axis.text.x = element_text(size=15,angle = 90))

#MEs
cmt=na.omit(match(mutburden[,1],rownames(clindata5)))
mutburden=mutburden[which(mutburden[,1] %in% rownames(clindata5)),]
memut=MEs[cmt,]
cln=clindata5[cmt,]
datExpr2=datExpr[,cmt]
prp=tempprop$Est.prop.weighted[cmt,]
#prp[prp==0]<-NA

mutcor=cor(mutburden$Freq,memut)
agecor=cor(cln$AGE,memut)
datcor=cor(mutburden$Freq,t(datExpr2))
propcor=cor(mutburden$Freq,prp,use='pairwise.complete.obs')

#Age Mutation Box Plot
data=data.frame((mutburden$Freq),cut(cln$AGE,c(0,25,35,45,55,65,75),c('20-29','30-39','40-49','50-59','60-69','70-79')))
#data=data.frame((mutburden$Freq),cut(cln$AGE,c(0,35,55,75),c('Young','Early-Aged','Late-Aged')))
names(data)<-c('Burden','Age')
data2=data%>%group_by(Age)%>%summarise(mn=mean(Burden),med=median(Burden))
m1=ggplot(data, aes(x = factor(Age),y = log10(Burden))) + geom_boxplot() + geom_point() + theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab('Age')+ylab('log Mutation Burden')
m1

#Genes that correlate with mutation burden
datcor=data.frame(gene=colnames(datcor),cor=datcor[1,])
df=datcor[order(datcor$cor),]
noquote(head(datcor$gene,50))
noquote(tail(datcor$gene,50))
datp<-p.adjust(WGCNA::corPvalueStudent(df$cor,307),method='fdr')
dfsig=df[datp<0.05,]

#Bulk RNAseq
coldata=data.frame(Mut=mutburden$Freq)
rownames(coldata)<-colnames(datExpr2)
dds <- DESeqDataSetFromMatrix(countData = datExpr2, colData = coldata, design = ~ Mut)

smallestGroupSize <- 100 # change from 3
keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize # change from 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)
#res<-res[res$log2FoldChange>0,]

resOrdered <- res[order(res$pvalue),]
resneg<-resOrdered[resOrdered$log2FoldChange<0,]
noquote(rownames(head(resOrdered,30)))

#Volcano
resshrink <- lfcShrink(dds, coef='Mut', res=res, type = 'normal')
yermyerm<-c('ABLIM1', 'MIB1', 'RNF19A', 'UBXN7', 'UBR3', 'UBE2H', 'UBE2K', 'UBE2E3', 'FBXL17', 'FBXO11', 'FBXO42', 'NPEPPS', 'VMP1', 'ULK1', 'RICTOR',
'SOS1', 'SOS2', 'MAP3K20', 'MAP4K4', 'MAPK8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND5', 'MT-CYB', 'VIPR1', 'IL7R', 'FCN3', 'HPGD',
'SEL1L', 'SMAD5', 'ANAPC1', 'UBE4A', 'USP25', 'USP33', 'RAD50', 'MTREX', 'IFI16', 'SAMHD1', 'CBL', 'APC', 'ETFB', 'COQ10A', 'MDH2', 'SCO2', 'IDH3B', 'COX4I1')
yermyerm=c(yermyerm,'ATF3', 'IER3', 'ING1', 'KDM6B', 'DUSP2', 'BIRC3', 'IRF1', 'IFITM1', 'NFKB2', 'MT-ND4L', 'DNAJC11', 'ATP13A2')#, 'ACAD9', 'ELAC2', 'SCO1')
yermyerm=c(yermyerm,rownames(resOrdered)[1:3],rownames(resneg)[1:3])
EnhancedVolcano(resshrink, lab=rownames(resshrink),selectLab = yermyerm, x = 'log2FoldChange', y = 'padj',pCutoff=0.05,FCcutoff = 2e-4,xlim=c(-0.001,0.001),drawConnectors = T,title=NULL,subtitle=NULL,caption=NULL)+theme_minimal()

#Compare with single-cell mutations
df=as.data.frame(resOrdered)
dfsig=df[df$pvalue<0.05,]
dffdr=df[df$padj<0.05,]

datal=read.table(paste0('/home/rd796/palmer_scratch/global/age-glmmTMB_global_newmut2_resultsTable.txt'),col.names=c('Gene','Beta0','Beta1','p_val','blank'),fill=T,skip=1)
datal=datal[!(is.na(datal$Beta1)|is.na(datal$p_val)|!is.na(datal$blank)),]
datal=mutate(datal,p_val_adj=p.adjust(datal$p_val,method='fdr',n=nrow(datal)))
datalsig=datal[datal$p_val<0.05,]
datalfdr=datal[datal$p_val_adj<0.05,]

intersect(datalsig$Gene,rownames(dfsig))
f1=length(which(datalfdr$Gene %in% rownames(dfsig))) #overlap
f2=length(which(datalfdr$Gene %in% rownames(df)))-f1 #in a not b
f3=length(which(rownames(dfsig) %in% datal$Gene))-f1 #in b not a
f4=length(which(rownames(df) %in% datal$Gene))-f1-f2-f3 #in neither
matrix <- matrix(c(f1, f2, f3, f4), nrow=2)
fisher.test(matrix, alternative="greater")

# Define the shared universe (genes tested in both experiments)
universe <- intersect(rownames(df), datal$Gene)
# Filter DE gene sets to include only genes in the shared universe
de_A <- intersect(datalfdr$Gene, universe)
de_B <- intersect(rownames(dfsig), universe)

f1=length(intersect(de_A, de_B))                       # In both A and B
f2=length(setdiff(de_A, de_B))                         # In A, not in B
f3=length(setdiff(de_B, de_A))                         # In B, not in A
f4=length(setdiff(universe, union(de_A, de_B)))        # In neither A nor B
matrix <- matrix(c(f1, f2, f3, f4), nrow=2)
fisher.test(matrix, alternative="greater")
chisq.test(matrix)

#SenMayo Mutation Scatter Plot
library(ggrepel)
data=data.frame(t(mutcor),t(agecor))
data=data %>% mutate(data,nm=substring(rownames(data),3))
names(data)<-c('MutCor','AgeCor','Name')
m2=ggplot(data,aes(x=AgeCor,y=MutCor)) + geom_point() + geom_label_repel(aes(label=ifelse(Name %in% c('2','15','19','20','28','31') ,Name,'')))+theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab('Age Correlation')+ylab('Mutation Burden Correlation')+
geom_hline(yintercept=0,linetype='dashed',color='red')+geom_vline(xintercept=0,linetype='dashed',color='red')
m2

#Which cell types have the most mutations
propcor=data.frame(colnames(propcor),propcor[1,])

#FOREST
library(forestplot)
library(DescTools)
library(stats)
TraitPvalue = corPvalueStudent(propcor, nrow(prp))
nms=colnames(propcor)
ll=matrix(0,3,ncol(prp))
for (x in 1:ncol(prp)) {
    ll[,x]=CorCI(propcor[x], nrow(prp), conf.level = 0.95, alternative = c("two.sided"))
}
fz=FisherZ(propcor)
f=order((fz),decreasing=TRUE)
lbl=cbind(nms[f],signif(TraitPvalue[f],1))
colnames(lbl)=c("Cell Type","p-value")
lbl=rbind(colnames(lbl),lbl)
ll=cbind(NA,ll)
f=c(1,f+1)
forestplot(col = fpColors(all.elements="black"),fn.ci_norm = fpDrawCircleCI,lwd.ci=2,lbl,t(ll)[f,],zero=0,is.summary = c(TRUE, rep(FALSE, 19)),xlab="Z-score",txt_gp = fpTxtGp(xlab=gpar(cex=1),ticks=gpar(cex=1)),boxsize=0.5,xticks=seq(from = -0.4, to = 0.4, by = 0.2))

#WGCNA HEATMAP
meow5=na.omit(cbind(memut,mutburden))
#moduleTraitCor = cor(meow5[,1:133], cbind(meow5[,134],meow5[,136]), use = "p")
cor1=cor(meow5[,1:133], meow5[,135], use = "p")
cor1 =cor1[order(as.numeric(substring(rownames(cor1),3))),]
cor1 =as.matrix(cor1)[2:133]
cor2=cor(MEs,clindata5$AGE,use="p")
cor2 =cor2[order(as.numeric(substring(rownames(cor2),3))),]
cor2 =as.matrix(cor2)[2:133]
cor3=cor2 #from plotscrna
moduleTraitCor = cbind(cor1,cor2,cor3)
moduleTraitPvalue = cbind(p.adjust(corPvalueStudent(cor1, nrow(meow5)),method="fdr"),p.adjust(corPvalueStudent(cor2, nrow(clindata5)),method="fdr"),p.adjust(corPvalueStudent(cor2, nrow(clindata5)),method="fdr"))
moduleTraitPvalue=cbind(moduleTraitPvalue,fetlist[,6])
moduleTraitPvalue=cbind(moduleTraitPvalue,matrix(1,132,1))
rownames(moduleTraitPvalue)=paste("ME",(1:(ncol(MEs)-1)))
cond=((moduleTraitPvalue[,1]<0.05|moduleTraitPvalue[,2]<0.05|moduleTraitPvalue[,4]<0.075)&holder[,1]!="null")
moduleTraitCor=moduleTraitCor[cond,]
moduleTraitPvalue=moduleTraitPvalue[cond,]
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
m=length(textMatrix)*3/5+1
n=length(textMatrix)*4/5
textMatrix[m:n]=paste("(",signif(moduleTraitPvalue[,4],2),")",sep="") #work on indices
m=length(textMatrix)*4/5+1
n=length(textMatrix)
textMatrix[m:n]=holder[cond,1] #work on indices
negs=moduleTraitCor<0
moduleTraitPvalue[,1:3]=moduleTraitPvalue[,1:3]
#moduleTraitPvalue=log10(abs(moduleTraitPvalue))
negs2=matrix(1,nrow(moduleTraitPvalue),ncol(moduleTraitPvalue))
negs2[negs]=-1
moduleTraitPvalue[,1:3]=moduleTraitPvalue[,1:3]*negs2[,1:3]

rampcols <- colorRampPalette(colors = c("white"), space="Lab")(200)
rampcols[91:100] <- colorRampPalette(colors = c("white","purple"), space="Lab")(10)
rampcols[101:110] <- colorRampPalette(colors = c("yellow","white"), space="Lab")(10)
xLabels=c('Mutation Burden','Age Correlation','sc_age','SenMayo','REACTOME Term')
m=length(textMatrix)*2/5+1
n=length(textMatrix)*3/5
rownames(moduleTraitPvalue)=str_replace(rownames(moduleTraitPvalue),"ME ","")
labeledHeatmap(Matrix = moduleTraitPvalue[,c(1,2,4,5)], xLabels = xLabels[-3], yLabels = rownames(moduleTraitPvalue), ySymbols = rownames(moduleTraitPvalue), 
colorLabels = FALSE, colors = rampcols, textMatrix = textMatrix[-m:-n], setStdMargins = TRUE, cex.text = 0.5, cex.lab.x=1,zlim = c(-1,1), 
main = paste("Module Attributes"))
