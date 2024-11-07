library('data.table')
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
options(bitmapType='cairo')

x2=data.table(fread("ltrcdata.txt",skip=31,header=FALSE,fill=TRUE)) #load ltrc data
x3=t(x2[c(-1:-9,-18:-43),])
colnames(x3)=x3[1,]
x3=x3[-1,]
x3=x3[x3[,1]=="disease state: Control",]
rownames(x3)=x3[,9]
x3=x3[,c(-1,-9)]
lagedata=as.data.frame(apply((x3[,-1:-7]),2,as.numeric))
rownames(lagedata)=rownames(x3)
lclindata=as.data.frame(x3[,1:7])
lclindata[,1]=as.numeric(substr(lclindata[,1],6,6))-1
lclindata[,2]=as.numeric(substr(lclindata[,2],6,7))
colnames(lclindata)=c('SEX','AGE','1','2','3','4','5')

annot = data.table(fread("ltrcannot.txt",skip=17, header=TRUE,fill=TRUE))
probes=colnames(lagedata)
probes2annot = match(probes, annot$ID)
colnames(lagedata)=annot$GENE_SYMBOL[probes2annot]

load(file = "evennewerstep1.RData")
lnames = load(file = "loggedmes.RData")
datExpr=log2(agedata3+1)
datTraits=clindata5
ids=read.csv('idgenename.csv')
probes = colnames(datExpr)
probes2annot = match(probes, ids$id)
colnames(datExpr) = ids$Description[probes2annot]
multiLabels=as.data.frame(moduleLabels)

probes2annot=unique(match(colnames(lagedata),colnames(datExpr)))
datExpr=datExpr[,na.omit(probes2annot)]
multiLabels=as.data.frame(multiLabels[na.omit(probes2annot),])
rownames(multiLabels)=colnames(datExpr)
probes2annot=unique(match(colnames(datExpr),colnames(lagedata)))
lagedata=lagedata[,na.omit(probes2annot)]
lagedata=lapply(lagedata, as.numeric)

multiColor=list(Female=t(multiLabels))
setLabels=c("GTEx","LTRC")
multiExpr = list(Female = list(data = as.data.frame(datExpr)), Male = list(data = as.data.frame(lagedata)))
system.time({mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 200, randomSeed = 1, quickCor = 0, verbose = 3)})

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][,-1], mp$preservation$observed[[ref]][[test]][,-1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][,-1], mp$preservation$Z[[ref]][[test]][,-1])
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)))
modColors = rownames(mp$preservation$observed[[ref]][[test]]) moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]
plotMods = modColors #%in% c('2','31','28','15','19','20')
text = modColors[plotMods]
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary")
par(mfrow = c(1,2)) 
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2) {
    min = min(plotData[, p], na.rm = TRUE)
    max = max(plotData[, p], na.rm = TRUE)
if (p==2) { if (min >-max/10) min =-max/10
ylim = c(min- 0.1 * (max-min), max + 0.1 * (max-min)) 
} else {ylim = c(max + 0.1 * (max-min), min- 0.1 * (max-min))}
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21, main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
if (p==2) { 
abline(h=0) 
abline(h=2, col = "blue", lty = 2) 
abline(h=10, col = "darkgreen", lty = 2) }
}
