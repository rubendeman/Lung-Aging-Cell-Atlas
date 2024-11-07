setwd('/home/rd796/project/ageproj')
library(WGCNA)
library(data.table)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
options(bitmapType='cairo')

x2=data.table(fread("ltrcdata.txt",skip=31, header=FALSE,fill=TRUE)) #load ltrc data
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
datExpr=log2(agedata3+1)
datTraits=clindata5
ids=read.csv('idgenename.csv')
probes = colnames(datExpr)
probes2annot = match(probes, ids$id)
colnames(datExpr) = ids$Description[probes2annot]

probes2annot=unique(match(colnames(lagedata),colnames(datExpr)))
datExpr=datExpr[,na.omit(probes2annot)]
probes2annot=unique(match(colnames(datExpr),colnames(lagedata)))
lagedata=lagedata[,na.omit(probes2annot)]
lagedata=lapply(lagedata, as.numeric)
multiExpr = vector(mode = "list", length = 2)
multiExpr[[1]] = list(data = as.data.frame(datExpr)) 
names(multiExpr[[1]]$data) = names(datExpr) 
rownames(multiExpr[[1]]$data) = rownames(datExpr) 
multiExpr[[2]] = list(data = as.data.frame((lagedata))) 
names(multiExpr[[2]]$data) = names(lagedata) 
rownames(multiExpr[[2]]$data) = rownames(lagedata) 
exprSize = checkSets(multiExpr)

#net = blockwiseModules(as.data.frame(multiExpr[[2]]), power = 6, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3)
net = blockwiseConsensusModules(multiExpr, power = 6, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3)
setLabels = c("GTEx", "LTRC")
plotEigengeneNetworks(net$multiMEs, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

concol=as.matrix(net$colors)
#rownames(concol)=substring(rownames(concol),6)#delete for multi
#conMEs=net$MEs
conMEs=as.data.frame(net$multiMEs[[1]])
nmes=ncol(conMEs)
lnames = load(file = "loggedmes.RData") 
gtMEs=MEs
gtcol=as.matrix(moduleLabels)
ids=read.csv('idgenename.csv')
probes = rownames(gtcol)
probes2annot = match(probes, ids$id)
tempnames=ids$Description[probes2annot]
probes2annot=unique(match(rownames(concol),tempnames))
gtcol=as.data.frame(gtcol[na.omit(probes2annot),])
rownames(gtcol)=tempnames[na.omit(probes2annot)]
probes2annot=unique(match(rownames(gtcol),rownames(concol)))#del
concol=as.data.frame(concol[na.omit(probes2annot),])#del
pTable = matrix(0, nrow = 133, ncol = nmes) 
for (fmod in 1:133){ 
    for (cmod in 1:nmes) { 
        #pTable[fmod, cmod] =sum(gtcol == as.numeric(substring(colnames(gtMEs)[fmod],3)) & concol == as.numeric(substring(colnames(conMEs[[1]])[cmod],3)))
        pTable[fmod, cmod] =sum((gtcol == fmod) & (concol == cmod))
    }
}
temp=rowSums(pTable)
rownames(pTable)=c(1:133)
colnames(pTable)=c(1:nmes)
pTable=pTable[temp>5,]
Heatmap(t(scale(t(pTable))),col=circlize::colorRamp2(c(0, 6), c("white", "red")),cluster_rows=FALSE,cluster_columns=FALSE)    
pTable=pTable[(rownames(pTable) %in% c('19','20','15','28','2','31')),]
#pTable=pTable[(rownames(pTable) %in% c('2','20','3','107','6','1')),]
Heatmap(t(scale(t(pTable))),col=circlize::colorRamp2(c(0, 6), c("white", "red")),cluster_rows=FALSE,cluster_columns=FALSE)

# Calculate the correlations
Traits=list()
Traits[[1]]=datTraits
Traits[[2]]=lclindata[,-3:-7]
moduleTraitCor = list()
moduleTraitPvalue = list()
for (set in 1:2){
  moduleTraitCor[[set]] = cor(as.data.frame(net$multiMEs[[set]]), Traits[[set]]$AGE, use = "p")
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], nrow(Traits[[set]]))
}
set=2
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), sep = "")
labeledHeatmap(Matrix = moduleTraitCor[[set]], xLabels = "AGE", yLabels = rownames(moduleTraitCor[[set]]), ySymbols = rownames(moduleTraitCor[[set]]),
colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, cex.lab=1, zlim = c(-1,1),
main = paste("Module--trait relationships"))
newset=as.data.frame(moduleTraitCor[[set]][(substring(rownames(moduleTraitCor[[set]]),8) %in% c('18','16','7','11','4','21')),])
newset=cbind(newset,gcor[match(c('2','31','28','15','19','20'),rownames(gcor)),])#Add other GTEx Cor
rownames(newset)=c('2','31','28','15','19','20')
textMatrix = signif(newset[,1:2], 2)
labeledHeatmap(Matrix = as.matrix(newset), xLabels = c("LTRC","GTEx"), yLabels = rownames(newset), ySymbols = rownames(newset),
colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, cex.lab=1, zlim = c(-1,1),
main = paste("Module--trait relationships"))

#consensus analysis
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])
textMatrix =  paste(signif(consensusCor, 2), "\n(", signif(consensusPvalue, 1), ")", sep = "")
labeledHeatmap(Matrix = consensusCor, xLabels = "AGE", yLabels = rownames(moduleTraitCor[[set]]), ySymbols = rownames(moduleTraitCor[[set]]),
colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, cex.lab=1, zlim = c(-1,1),
main = paste("Module--trait relationships"))

datTraits=as.data.frame(lclindata[,-3:-7])
datExpr=as.data.frame(lagedata)
MEs=as.data.frame(net$multiMEs[[2]])
names(MEs) = substring(names(MEs), 6)
moduleLabels=net$colors
