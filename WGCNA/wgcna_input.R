setwd('/home/rd796/project/ageproj')
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
options(bitmapType='cairo')

agedata = read.csv("mat2.csv",stringsAsFactors = FALSE)
agedata2=agedata[,-c(1)]
agedata3 = as.data.frame(t(agedata2))
colnames(agedata3)=agedata$ID_REF
gsg = goodSamplesGenes(agedata3, verbose = 3)
sampleTree = hclust(dist(agedata3), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
clust = cutreeStatic(sampleTree, cutHeight = 120)
table(clust)
keepSamples = (clust==1)
agedata4 = agedata3[keepSamples, ]

clindata = read.csv("clin2.csv",stringsAsFactors = FALSE)
clindata2=clindata[,-c(1)]
clindata3 = as.data.frame(t(clindata2))
colnames(clindata3)=clindata$ID_REF
traitRows = match(rownames(agedata4), rownames(clindata3))
clindata4 = clindata3[traitRows,]

sampleTree2 = hclust(dist(agedata4), method = "average")
traitColors = numbers2colors(clindata4, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(clindata4), main = "Sample dendrogram and trait heatmap")
save(agedata4, clindata4, file = "newstep1.RData")
checker = match(rownames(agedata4), rownames(clindata4))
checker
