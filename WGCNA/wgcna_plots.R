setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(WGCNA)

load("loggedmes.RData") #load module eigengenes
load("evennewerstep1.RData") #load bulk
load("geneInfo0.RData") #load module genes
datExpr=log2(agedata3+1)
datTraits=clindata5
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

#HOW MANY COR AGE
gcor=cor(MEs,clindata5$AGE,use="p")
gcor =gcor[order(as.numeric(substring(rownames(gcor),3))),]
gcor =as.matrix(gcor[2:133])
gpv=p.adjust(corPvalueStudent(gcor, nrow(clindata5)),method="fdr")
gpv=as.data.frame(gpv)
gcor=as.data.frame(gcor)
#rownames(cor2)=paste0("ME",(1:(ncol(MEs)-1)))
rownames(gcor)=paste((1:(ncol(MEs)-1)))
rownames(gpv)=paste((1:(ncol(MEs)-1)))
agecorgenes=as.numeric(rownames(gpv)[gpv<0.05])

poscorgenes=as.numeric(rownames(gpv)[gpv<0.05&gcor>0])
negcorgenes=as.numeric(rownames(gpv)[gpv<0.05&gcor<0])

#Plot module size vs age cor
modsizes=unlist(sapply(1:(ncol(MEs)-1),function(a)length(which(as.vector(geneInfo0$moduleColor)==a))))
plot(as.vector(modsizes), as.vector(t(-log10(gpv))), col = 1,  pch = 21, main = 'Module Size vs Significance', cex = abs(as.vector(t(-log10(gpv))))*1, ylab = '-log(p-value)', xlab = "Module Size", log = "x", cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(as.vector(modsizes)[gpv<0.05], as.vector(t(-log10(gpv)))[gpv<0.05], rownames(gpv)[gpv<0.05], cex = 1, offs = 0.06)

#FIGURE 1
library(ggrepel)
df=as.data.frame(cbind(as.vector(log10(modsizes)), as.vector(t(-log10(gpv))),as.vector(t(-log10(gpv)*sign(gcor)))))
colnames(df)=c('Genes','Significance','Correlation')
g3<-ggplot(df, aes(Genes,Significance,label=rownames(df))) + geom_point(aes(color = Correlation,size = Significance), shape = 16)  + labs(x="log Module Size",y="-log10(p-value)",color='Significance',size='-log10(p-value)') + 
scale_x_continuous(expand=c(0,0),limits=c(0,6))+scale_y_continuous(expand=c(0,0),limits=c(0,12)) + geom_label_repel(aes(label=ifelse(Significance>3.9,rownames(df),'')),label.size=NA,fill=NA) + 
scale_color_gradient2(low="purple", mid="gray",high="yellow", space ="Lab") + theme_minimal() + theme(legend.position = c(0.85,0.5),axis.line=element_line(size=1),panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
g3

#Bulk age cor vs scRNA age cor
plot(cor2,cor3,pch = 19, main = 'Bulk age correlation VS scRNAseq age correlation', ylab = 'scRNAseq age correlation', xlab = 'Bulk age correlation',cex.main =0.9)
text(cor2[c(2,15,19,20,28,31)],cor3[c(2,15,19,20,28,31)],names(cor2)[c(2,15,19,20,28,31)],cex=1,pos=3)

#Age cor vs mutcor
par(mar = c(4,4,4,2))
names(cor2)=1:132
plot(cor2[1:132],cor1[1:132],pch = 19, main = 'Age correlation VS mutation burden correlation', ylab = 'Global mutation correlation', xlab = 'Age correlation',cex.main =0.9)
text(cor2[c(2,15,19,20,28,31)],cor1[c(2,15,19,20,28,31)],names(cor2)[c(2,15,19,20,28,31)],cex=1,pos=3)

#All muts vs mod muts
par(mar = c(4,4,4,2))
plot(modcorlist[2:133,1],cor1[1:132],xlim=c(-0.2,0.2), ylim=c(-0.2,0.2),pch = 19, main = 'Intra-module mutation correlation VS global mutation correlation', ylab = 'Global mutation correlation', xlab = 'Intra-module mutation correlation',cex.main =0.9)
text(modcorlist[c(3,20,21,16,29,32),1],cor1[c(2,15,19,20,28,31)],names(cor2)[c(2,15,19,20,28,31)],cex=1,pos=3)

#3D
library("scatterplot3d")
s=scatterplot3d(cor2[1:132],cor1[1:132],modcorlist[2:133,1],pch=19,angle=70,main = 'Age correlation VS mutation burden correlation',xlab='Age correlation',ylab='Global mutation burden',zlab='Module mutation burden',cex.main =0.9)
#ind=(plothub[,4]<0.1)
text(s$xyz.convert(cor2[c(2,15,19,20,28,31)],cor1[c(2,15,19,20,28,31)],modcorlist[c(3,16,20,21,29,32),1]),names(cor2)[c(2,15,19,20,28,31)],cex=1,pos=3)

#GO Enrichment Heatmap Modules
library('gprofiler2')
lnames=load('geneInfo0.RData')
holder=matrix(0,18,1)
counter=-2
for (i in c(2,15,19,20,28,31)){
  counter=counter+3
temp=geneInfo0$ensg[geneInfo0$moduleColor==i]
k=gost(query=temp, organism = "hsapiens")
if(is.null(k)){
term="null"
}else{
term=k$result$term_id[k$result$source=="REAC"]
}
holder[counter:(counter+2),1]=term[1:3]
}
temp1=temp=geneInfo0$ensg[geneInfo0$moduleColor==2]
temp2=temp=geneInfo0$ensg[geneInfo0$moduleColor==15]
temp3=temp=geneInfo0$ensg[geneInfo0$moduleColor==19]
temp4=temp=geneInfo0$ensg[geneInfo0$moduleColor==20]
temp5=temp=geneInfo0$ensg[geneInfo0$moduleColor==28]
temp6=temp=geneInfo0$ensg[geneInfo0$moduleColor==31]
k=gost(query=list(temp1,temp2,temp3,temp4,temp5,temp6), multi_query=TRUE, organism = "hsapiens")
result=k$result
result2=result[result$term_id %in% holder,]
result3=-log10(as.data.frame(result2$p_values))
rownames(result3)=c('2','15','19','20','28','31')
colnames(result3)=result2$term_name
#plot
ht=Heatmap(t(result3), name = "-log10(p)", col= colorRampPalette(colors = c("grey","blue"), space="Lab")(100),
           cluster_columns = FALSE, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = TRUE,
           show_row_dend = FALSE, row_names_gp = gpar(fontsize = 10),heatmap_legend_param = list(legend_direction = "horizontal"))
draw(ht, padding = unit(c(8, 2, 2, 12), "mm"),heatmap_legend_side = "bottom")
#or invert
ht=Heatmap((result3), name = "-log10(p)", col= colorRampPalette(colors = c("grey","blue"), space="Lab")(100),
           cluster_columns = T, show_column_dend = FALSE, column_title_rot = 45, cluster_rows = F,
           show_row_dend = FALSE, row_names_gp = gpar(fontsize = 10),heatmap_legend_param = list(legend_direction = "vertical"))
draw(ht, padding = unit(c(30, 2, 2, 12), "mm"),heatmap_legend_side = "right")

#all age mods
library(gprofiler2)
glist=list()
for (i in agecorgenes){ #negcorgenes
  glist[[i]]=geneInfo0$ensg[geneInfo0$moduleColor==i]
}
k1=gost(query=glist, multi_query=FALSE, organism = "hsapiens")
kres=k1$result[k1$result$source=='REAC',]
kres=kres[!duplicated(kres$query),]
kmat=kres[,c(1,3,11)]
kmat[,1]=substring(kmat[,1],7)
kmat=kmat %>% mutate(logp=-log10(gpv[kmat$query,]),sign=sign(gcor[kmat$query,]))
colnames(kmat)=c('Module','p_value','term_name','logp','sign')
kmat=kmat[order(kmat$logp,decreasing=FALSE),]
kmat=kmat %>% mutate(place = if_else(row_number() >21, 1, 0)) 
kmat=kmat %>% mutate(color = if_else(row_number() >21, 'white', 'black')) 
kmat[19,3]='TCA cycle and respiratory electron transport'
kmat[6,3]='Immunoregulatory interactions'
kmat=kmat[-22,]
g4=ggplot(kmat, aes(x = logp, y = Module, fill=logp*sign(sign))) + geom_col()+scale_fill_gradient2(low="purple", mid="gray", high="yellow",midpoint=0)+
    geom_text(aes(label = term_name),hjust = kmat$place,colour=kmat$color) + scale_y_discrete(limits=kmat$Module) +
    theme_minimal()+theme(axis.line.x=element_line(size=1),legend.position=c(0.85,0.25),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(x='-log10(p-value)',fill='-log10(p-value)')+NoLegend()
    g3+g4

#FET Testing
sm=read.table('listsenmayo.txt') 
colnames(sm)='gene_ID'
gen=read.csv('GenAge.csv')
colnames(gen)[2]='gene_ID'
cell=read.table('CellAge.csv',sep=';',header=TRUE)
colnames(cell)[2]='gene_ID'
fridman=read.csv('FRIDMAN.csv')
purcell=read.csv('SEGURA.csv')
segura=read.csv('PURCELL.csv')
datalist=list(cell,gen,purcell,fridman,segura,sm)
datacontainer=matrix(0,ncol(MEs)-1,6)
fetlist=matrix(0,ncol(MEs)-1,6)
for (i in 1:(ncol(MEs)-1)){
temp=geneInfo0$geneSymbol[geneInfo0$moduleColor==i]
for (j in 1:6){
  slist=datalist[[j]]$gene_ID
cnt=temp %in% slist
datacontainer[i,j]=sum(cnt)
#fetobj=matrix(c(sum(cnt),length(slist),length(temp),54073),2,2)
fetobj=matrix(c(sum(cnt),length(slist)-sum(cnt),length(temp)-sum(cnt),54073-(length(slist)-sum(cnt))-(length(temp)-sum(cnt))-sum(cnt)),2,2)
fet=fisher.test(fetobj,alternative="greater")
fetlist[i,j]=fet$p.value
}
}
rownames(fetlist)=as.character(1:(ncol(MEs)-1))
cond=(fetlist[,1]<0.05|fetlist[,2]<0.05|fetlist[,3]<0.05|fetlist[,4]<0.05|fetlist[,5]<0.05|fetlist[,6]<0.05)
fetlist=fetlist[cond,]
datacontainer=datacontainer[cond,]

#SUPPLEMENT: SENESCENCE LIST HEATMAP
rampcols <- colorRampPalette(colors = c("white"), space="Lab")(100)
rampcols[1:5] <- colorRampPalette(colors = c("purple","white"), space="Lab")(5)
labeledHeatmap(Matrix = fetlist, ylab="Module", xLabels = c('CellAge','GenAge','PURCELL','FRIDMAN','SEGURA','SenMayo'), yLabels=rownames(fetlist),
   colorLabels = FALSE, colors = rampcols, textMatrix = datacontainer, setStdMargins = FALSE, cex.text = 0.5,legendLabel = "p-value", 
   zlim = c(0,1), main = NULL)

#SUPPLEMENT: SENESCENCE LIST UPSET
library('UpSetR')
slist=list()
for (j in 1:6){slist[[j]]=datalist[[j]]$gene_ID}
names(slist) <- c('CellAge','GenAge','PURCELL','FRIDMAN','SEGURA','SenMayo')
upset(fromList(slist), nsets=6,nintersects=17,order.by=c('freq','degree'),decreasing=c(TRUE,FALSE),text.scale=2)
