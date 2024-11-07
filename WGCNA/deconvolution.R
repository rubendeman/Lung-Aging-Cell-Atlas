setwd('/home/rd796/project/ageproj')
library(MuSiC)
library(cowplot)
library(Biobase)

load('imm12_10.RData') #load integrated

#PERFORM FILTERING AS IN VAL_SCRNA
decon.ref<-immune.combined
#decon.ref<-subset(immune.combined,subset=Manuscript_Identity=='IPF Cell Atlas')

#Downsample
all.cellTypes <- unique(decon.ref@meta.data$predicted.id)
all.cells.remove=list()

#Loop
for(i in 1:length(all.cellTypes)){
    temp.meta <- decon.ref@meta.data %>% filter(predicted.id==all.cellTypes[i])
    nCells <- nrow(temp.meta)
    if(nCells > 500){
        ind=(1:nrow(decon.ref@meta.data))[decon.ref@meta.data$predicted.id==all.cellTypes[i]]
        cells.remove <- sample(ind, size=(nCells-500), replace=FALSE)
    }
all.cells.remove[[i]]<-cells.remove
}
decon.ref<- decon.ref[,-unlist(all.cells.remove)]

all(rownames(decon.ref@meta.data) == colnames(decon.ref@assays$RNA@counts))

phenoData <- new("AnnotatedDataFrame", data=decon.ref@meta.data)

decon.ref <- ExpressionSet(assayData=as.matrix(decon.ref@assays$RNA@counts), phenoData=phenoData) #use raw count data, not normalized (recommended in vignette)

################
# Work on bulk #
################
load(file = "evennewerstep1.RData")
datExpr=agedata3
datTraits=clindata5;
ids=read.csv('idgenename.csv')
probes = colnames(datExpr)
probes2annot = match(probes, ids$id)
colnames(datExpr) = ids$Description[probes2annot]

bsm.counts <- t(datExpr)
bsm.design <- datTraits
colnames(bsm.counts) <- rownames(bsm.design)

dup.genes <- unique(rownames(bsm.counts)[duplicated(rownames(bsm.counts))])

rows.remove <- which(rownames(bsm.counts) %in% dup.genes)
bsm.counts <- bsm.counts[-c(rows.remove),]

### clean it up and make it a matrix
bsm.counts <-  data.matrix(bsm.counts)

#### okay now make it an expressionSet
bsm.phenoData <- new("AnnotatedDataFrame", data=bsm.design)
bsm.bulk.ExpSet <- ExpressionSet(assayData=bsm.counts, phenoData=bsm.phenoData)

est.prop <- music_prop(bulk.eset=bsm.bulk.ExpSet,
    sc.eset=decon.ref, 
    clusters="predicted.id",
    samples="orig.ident")
