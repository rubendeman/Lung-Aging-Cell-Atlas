setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(WGCNA)
#options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
#options(bitmapType='cairo')

load(file = "loggedmes.RData") #load module eigengenes
load(file = "evennewerstep1.RData") #load bulk
datExpr=log2(agedata3+1)
datTraits=clindata5
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)

#FOR MODULES
TraitCor = cor(MEs, datTraits, use = "p")
TraitPvalue = corPvalueStudent(TraitCor, nSamples)
moduleTraitPvalue=data.frame()
moduleTraitCor=data.frame()
library(ppcor)
  for(i in 1:(length(table(moduleLabels)))){
    x=na.omit(data.frame(cbind(MEs[,i],datTraits$AGE,datTraits$SEX)))
    x=pcor.test(x[,1],x[,2],x[,3], method="spearman")
    moduleTraitPvalue[i,1]=x$p.value
    moduleTraitCor[i,1]=x$estimate
    }
moduleTraitPvalue=p.adjust(as.matrix(moduleTraitPvalue), method="fdr")
rownames(moduleTraitPvalue)=names(MEs)
rownames(moduleTraitCor)=names(MEs)

trait.est=NULL
trait.pval=NULL
  for(i in 1:(length(table(moduleLabels)))){
    x=na.omit(data.frame(cbind(gene=MEs[,i],
                               SEX=datTraits$SEX,
                               AGE=datTraits$AGE)))
    lm.model=lm(gene~AGE+SEX+AGE:SEX, data=x)
  trait.pval=rbind(trait.pval,summary(lm.model)$coefficients[,4])
  trait.est=rbind(trait.est,summary(lm.model)$coefficients[,1])
    }
trait.pval=p.adjust(trait.pval, method='fdr')
rownames(trait.est)=names(MEs)
rownames(trait.pval)=names(MEs)

corrchart=cbind(TraitCor, TraitPvalue, moduleTraitCor, moduleTraitPvalue, trait.est, trait.pval)

#write.csv(corrchart, file = "lmoutput.csv")

#FOR GENES
nSamples = nrow(datExpr)
weight = as.data.frame(datTraits$AGE)
names(weight) = "age"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

ids=read.csv('idgenename.csv')
probes = names(datExpr)
probes2annot = match(probes, ids$id)
geneInfo0 = data.frame(substanceBXH = probes, ensg = ids$Name.1[probes2annot], geneSymbol = ids$Description[probes2annot],
moduleColor = moduleLabels, geneTraitSignificance, GSPvalue)
modOrder = order(-abs(cor(MEs, weight, use = "p")))

for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]])
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
#Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.age))
geneInfo = geneInfo0[geneOrder, ]

#write.csv(geneInfo, file = "loggeneInfo.csv")
