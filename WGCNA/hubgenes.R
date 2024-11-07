setwd('/home/rd796/project/ageproj')
library(WGCNA)

load('geneInfo0.RData') #load module genes

#Load Senescence Lists
sm=read.table('listsenmayo.txt') 
colnames(sm)='gene_ID'
gen=read.csv('GenAge.csv')
colnames(gen)[2]='gene_ID'
cell=read.table('CellAge.csv',sep='',header=TRUE)
colnames(cell)[2]='gene_ID'
fridman=read.csv('FRIDMAN.csv')
purcell=read.csv('SEGURA.csv')
segura=read.csv('PURCELL.csv')
allsen=c(cell$gene_ID,gen$gene_ID,purcell$gene_ID,fridman$gene_ID,segura$gene_ID,sm$gene_ID)

i=1
gs=geneInfo0$GS.age[geneInfo0$moduleColor==i]
nm=geneInfo0$geneSymbol[geneInfo0$moduleColor==i]
whichcol=as.numeric(colnames(geneInfo0) %>% `==`(paste("MM.",i,sep="")) %>% which())
mm=geneInfo0[geneInfo0$moduleColor==i,whichcol]
hub=cbind(nm,gs,mm)
hub=cbind(hub,abs(gs*mm))
hub=cbind(hub,hub[,1] %in% allsen)
hub=hub[order(hub[,4],decreasing=TRUE),]
View(hub)

library('gprofiler2')
k=gost(query=head(hub[,1],30), organism = "hsapiens")
View(k$result)
