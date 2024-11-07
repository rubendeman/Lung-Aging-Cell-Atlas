setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(ggplot2)
library(dplyr)
library(WGCNA)
library(FedData)
load('evennewerstep1.RData') #load bulk
load('loggedmes.RData') #load module eigengenes
load('taball.RData') #load bulk mutation data
load('prop12_15.RData') #load music props
head(tab_all,3)
tbl=tab_all[tab_all$tissue == "Lung",]
mutburden=as.data.frame(table(tbl$gtexIds))

comp<-list(epi=c('AT1','AT2','Basal','Ciliated','Club','Goblet'),endo=c('Lymphatic','Peribronchial','Aerocyte','gCap','Arterial','Venous'),mes=c('Adventitial Fibroblast','Alv. Fibroblast','Myofibroblast','SMC','Pericyte'),imm=c('Monocyte','Macrophage','Alv. Macrophage'),lymph=c('B','T','Mast','DC','NK'))
cell.types<-unlist(comp)

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
prp=tempprop$Est.prop.weighted[cmt,]
#prp[prp==0]<-NA

mutcor=cor(mutburden$Freq,memut)
agecor=cor(cln$AGE,memut)
propcor=cor(mutburden$Freq,prp,use='pairwise.complete.obs')

#Age Mutation Box Plot
data=data.frame((mutburden$Freq),cut(cln$AGE,c(0,25,35,45,55,65,75),c('20-29','30-39','40-49','50-59','60-69','70-79')))
names(data)<-c('Burden','Age')
m1=ggplot(data, aes(x = factor(Age),y = log10(Burden))) + geom_boxplot() + geom_point() + theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab('Age')+ylab('log Mutation Burden')
m1

#SenMayo Mutation Scatter Plot
library(ggrepel)
data=data.frame(t(mutcor),t(agecor))
data=data %>% mutate(data,nm=substring(rownames(data),3))
names(data)<-c('MutCor','AgeCor','Name')
m2=ggplot(data,aes(x=AgeCor,y=MutCor)) + geom_point() + geom_label_repel(aes(label=ifelse(Name %in% c('2','15','19','20','28','31') ,Name,'')))+theme_minimal()+theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlab('Age Correlation')+ylab('Mutation Burden Correlation')+
geom_hline(yintercept=0,linetype='dashed',color='red')+geom_vline(xintercept=0,linetype='dashed',color='red')
m2

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
