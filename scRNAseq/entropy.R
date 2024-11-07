setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(philentropy)

# ENTROPY FUNCTION
sc.ent.KL.count <- function(X){
  X <- unlist(as.vector(X))
  null.x <- rep(1,length(X))
  #null.x <- null.x/sum(null.x)
  # check
  #sum(X$value) # should equal 1
  #sum(null.x) # should equal 1
  KL.count <- philentropy::KL(x=rbind(X,null.x),# null probability distribution
                          est.prob = 'empirical', # these are already probability distributions that add to 1, so we will not treat as counts
                          unit = 'log' # use nats as output unit
                          )
  return(KL.count)
}
rscale1=function(x){(x-min(x))/(max(x)-min(x))}
rscale2=function(x){1-(x-min(x))/(max(x)-min(x))}

# Define the celltypes ####
celltypes <- unique(immune.combined@meta.data$predicted.id)
#celltypes <- c('AT1','AT2B','AT2S','gCap','Aerocyte')
celltypes <- c('AT1','AT2','gCap','Aerocyte')

### Downsample Method
Idents(immune.combined)<-paste(immune.combined$predicted.id)
seurat.object<-subset(immune.combined,downsample=10000)

# Wrapper function to calculate entropy for a given cell type
calculateEntropy <- function(celltype){
  print(paste("Working on", celltype))
  tmp <- subset(seurat.object, subset=predicted.id==celltype)
  expr <- tmp@assays$RNA@counts

  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])

  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)

  nsample <- min(table(tmp@meta.data$agebin)[c("Young", "Aged")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$agebin == "Aged")], nsample)
young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$agebin == "Young")], nsample)
    ds_expr_r <- ds_expr[, c(young_r, old_r)]
    genes <- rownames(ds_expr_r)

    calcEuclDist <- function(matr, young, old){
    tmp <- data.matrix(sqrt(matr[genes, young]))
    d_young <- apply(tmp,2,function(x){sc.ent.KL.count(x)})

    tmp <- data.matrix(sqrt(matr[genes, old]))
    d_old <- apply(tmp,2,function(x){sc.ent.KL.count(x)})

    list(Young = d_young, Aged = d_old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, young = young_r, old = old_r)
  ds
}

# Run for all celltypes ####
res <- lapply(celltypes, function(x) calculateEntropy(x))
names(res) <- celltypes

res_original <- res

# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "fdr")

################
# Generate Fig #
################
names(res)[names(res) %in% c('Alv. Fibroblast')]<-'Alv Fibroblast'
names(res)[names(res) %in% c('Alv. Macrophage')]<-'Alv Macrophage'

resm=unlist(res)
temp=data.frame(nm=names(resm),val=resm)
temp=cbind(data.frame(do.call('rbind', strsplit(as.character(temp$nm), '.', fixed=TRUE))),temp)
names(temp)=c('Type','Age','Ident','Barcode','Entropy')

#################################
# ADD METADATA (versus Entropy) #
#################################
mdata=data.frame(sid=seurat.object$orig.ident,n_age=seurat.object$age)
t2=cbind(temp,mdata[match(temp$Ident,rownames(mdata)),])

################
# Load results #
################
load('/home/rd796/scratch_public/Backup/Ruben/Entropy/entropyres_sqrt.RData')
load('fill_df.RData')

t2$Type[t2$Type=='AT2B']='SPCʰⁱᵍʰ AT2'
t2$Type[t2$Type=='AT2S']='SPCˡᵒʷ AT2'

#colors
mycolors=unlist(fill_df$color)
names(mycolors)=fill_df$predicted.id
mycolors=c(mycolors,'SPCʰⁱᵍʰ AT2'="5BB300",'SPCˡᵒʷ AT2'="5BB300")
names(mycolors)<-sub("\\.","",names(mycolors))
colfun <- scale_colour_manual(name='Cell Type',values=mycolors)

f=as.data.frame(sapply(unique(t2$Type),function(x){tryCatch({t.test(t2$Entropy[t2$Type==x&t2$Age=='Young'],t2$Entropy[t2$Type==x&t2$Age=='Aged'])},error=function(e){})}))
f['estimate',]
f['p.value',]

ggplot(t2, aes(y = rscale2(Entropy), x = factor(Type,levels=names(f)[order(unlist(f['statistic',]),decreasing=T)]), fill=factor(Age,levels=c('Young','Aged')))) + geom_boxplot() + geom_point(position=position_jitterdodge(jitter.width=0.1))+
xlab("Type")+ylab("Transcriptional Entropy (Scaled)")+scale_fill_manual(values=c('#91C4F2','#253C78'))+scale_y_continuous(limits=c(0,1.1))+
guides(fill=guide_legend(title=NULL))+theme_classic()+theme(axis.text.x = element_text(size=12,angle=90)) + stat_compare_means(label='p.signif')

ggplot(t2[t2$Type %in% c('gCap','SPCʰⁱᵍʰ AT2','Aerocyte','AT1','SPCˡᵒʷ AT2','Macrophage'),], aes(y = rscale2(Entropy), x = factor(Type,levels=names(f)[order(unlist(f['statistic',]),decreasing=T)]), fill=factor(Age,levels=c('Young','Aged')))) + geom_boxplot() + geom_point(position=position_jitterdodge(jitter.width=0.1))+
xlab("Type")+ylab("Transcriptional Entropy (Scaled)")+scale_fill_manual(values=c('#91C4F2','#253C78'))+scale_y_continuous(limits=c(0,1.1))+
guides(fill=guide_legend(title=NULL))+theme_classic()+theme(axis.text.x = element_text(size=12,angle=90)) + theme_classic() + theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+stat_compare_means(label='p.signif')

f=as.data.frame(sapply(unique(t2$Type),function(x){cor.test(t2$n_age[t2$Type==x],t2$Entropy[t2$Type==x],use='pairwise.complete.obs',alternative='less')}))
f['estimate',]
f['p.value',]

ggplot(t2, aes(y = Entropy, x = n_age, color = Type)) + colfun + geom_point() + theme_classic() + geom_smooth(method='lm',se=F) +
scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,6))

###################
# GENE EXPRESSION #
###################
ct=c('gCap','Arterial','AT1','AT2')
auto=c(gCap=list(data.frame(Gene=c('IL7R','FCN3','EDN1','APLNR','TEK','CA4','EDNRB','FENDRR','VIPR1','APLNR','APLN','HPGD'))),Arterial=list(data.frame(Gene=c('DKK2','GJA5','ACKR1'))),
AT1=list(data.frame(Gene=c('HOPX','AGER','RTKN2'))),AT2=list(data.frame(Gene=c('SFTPC','LAMP3','SLC34A2','SFTPB','SFTPA1'))))

gene=t(as.data.frame(seurat.object@assays$RNA@data[unlist(auto),]))
t2=cbind(t2,gene[match(temp$Ident,rownames(gene)),])

ggplot(t2[t2$Type %in% c('SPCʰⁱᵍʰ AT2','SPCˡᵒʷ AT2'),], aes(y = sftpc, x = rscale2(Entropy), color = Type)) + geom_point() + theme_classic() +
xlab('Entropy (Scaled)')+scale_color_manual(values=c('#5BB300','red'))+scale_x_continuous(limits=c(0,1.2))

#Check correlations without plotting
t2$Type[t2$Type %in% c('SPCʰⁱᵍʰ AT2','SPCˡᵒʷ AT2')]='AT2'
auto2=lapply(ct,function(y){cor(t2$Entropy[t2$Type==y],t2[t2$Type==y,unlist(auto[[y]])])})

#########################
# PSEUDOBULK PER SAMPLE #
#########################
t3<-t2 %>% group_by(sid, Type, Age) %>% select_if(is.numeric) %>% summarise_all(mean)

f=as.data.frame(sapply(unique(t3$Type),function(x){tryCatch({t.test(t3$Entropy[t3$Type==x&t3$Age=='Young'],t3$Entropy[t3$Type==x&t3$Age=='Aged'])},error=function(e){})}))
f['estimate',]
f['p.value',]

ggplot(t3, aes(y = Entropy, x = factor(Type,levels=names(f)[order(unlist(f['statistic',]),decreasing=T)]), fill=factor(Age,levels=c('Young','Aged')))) + geom_boxplot() + geom_point(position=position_jitterdodge(jitter.width=0.1))+
xlab("Type")+scale_fill_manual(values=c('#91C4F2','#253C78'))+scale_y_continuous(limits=c(0,5))+
guides(fill=guide_legend(title=NULL))+theme_classic()+theme(axis.text.x = element_text(size=12,angle=90)) + stat_compare_means(label='p.signif')

f=as.data.frame(sapply(unique(t3$Type),function(x){cor.test(t3$n_age[t3$Type==x],t3$Entropy[t3$Type==x],use='pairwise.complete.obs')}))
f['estimate',]
f['p.value',]

ggplot(t3, aes(y = Entropy, x = n_age, color = Type)) + colfun + geom_point() + geom_smooth(method='lm',se=F) + theme_classic() +
scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,5))

#####################
# SOMATIC MUTATIONS #
#####################
t3$Type[t3$Type %in% c('Alv Fibroblast','Adventitial Fibroblast')]<-'Fibroblast'

#load('cellmutmat.RData')
load('cellmutmat2.RData')

names(gfg_data)[names(gfg_data)=='AT2B']='SPCʰⁱᵍʰ AT2'
names(gfg_data)[names(gfg_data)=='AT2S']='SPCˡᵒʷ AT2'

sid=rownames(gfg_data)
gfg_data=gfg_data[,colSums(gfg_data,na.rm=T)!=0]
gfg_bin=apply(gfg_data,2,function(x){cut(x,c(0,median(x,na.rm=T),100),c('Low Mutation','High Mutation'))}) #optional, makes mutations binary
gfg_data<-pivot_longer(data.frame(sid,gfg_data),!sid)
gfg_bin<-pivot_longer(data.frame(sid,gfg_bin),!sid)
gfg_data$name=sub("\\.","",sub("\\."," ",gfg_data$name))
gfg_bin$name=sub("\\.","",sub("\\."," ",gfg_bin$name))
t4=cbind(t3,rateMut=gfg_data$value[match(paste0(t3$Type,t3$sid),paste0(gfg_data$name,gfg_data$sid))],binMut=gfg_bin$value[match(paste0(t3$Type,t3$sid),paste0(gfg_bin$name,gfg_bin$sid))])

library(broom)
f=sapply(unique(t4$Type),function(x){tryCatch({t.test(t4$Entropy[t4$Type==x&t4$binMut=='Low Mutation'],t4$Entropy[t4$Type==x&t4$binMut=='High Mutation'])},error=function(e){})})
f = f[-which(sapply(f, is.null))]
g=as.data.frame(lapply(f,function(x){unlist(tidy(x))}))
names(g)[names(g)=='Alv.Macrophage']='Alv Macrophage'
g['estimate',]
g['p.value',]

t4$Type=factor(t4$Type,levels=names(g)[order(unlist(g['p.value',]))])
ggplot(na.omit(t4), aes(y = Entropy, x = Type, fill=binMut)) + geom_boxplot() + geom_point(position=position_jitterdodge(jitter.width=0.1))+
xlab("Type")+scale_fill_manual(values=c('#A20021','#F52F57'))+scale_y_continuous(limits=c(0,5))+
guides(fill=guide_legend(title=NULL))+theme_classic()+theme(axis.text.x = element_text(size=12,angle=90)) + stat_compare_means(label='p.signif')

ggplot(na.omit(t4), aes(x = n_age, y = -log10(rateMut), color=Type)) + colfun + geom_point() + geom_smooth(method='lm',se=F) + 
xlab("Age") + ylab("Mutation Burden") + theme_classic() + scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,6))

##########
# GLOBAL #
##########
t5 <- t2 %>% group_by(sid) %>% select_if(is.numeric) %>% summarise(Entropy_mean=mean(Entropy),Entropy_sd=sd(Entropy))

#load('mutmat.RData')
load('mutmat2.RData')
mutmat=gfg_sample[gfg_sample$mut>0&gfg_sample$mut<3000,]
mutmat=mutmat %>% mutate(logMut=-log10(b),binMut=cut(b,c(0,0.00003,1),c('Low Mutation Burden','High Mutation Burden')))

temp=cbind(mutmat,t5[match(rownames(mutmat),t5$sid),2:3])
colnames(temp)<-c('Age','Calls','nMut','rateMut','logMut','binMut','Entropy','sdEntropy')

f=t.test(temp$Entropy[temp$binMut=='Low Mutation Burden'],temp$Entropy[temp$binMut=='High Mutation Burden'])
f['estimate']
f['p.value']

ggplot(temp, aes(x = binMut, y = Entropy, fill = binMut)) + geom_boxplot() + geom_point(position=position_jitterdodge(jitter.width=0.1)) +
xlab("Type")+ylab("Organization")+scale_fill_manual(values=c('#F52F57','#A20021'))+scale_y_continuous(limits=c(0,5))+
guides(fill=guide_legend(title=NULL))+theme_classic()+theme(axis.text.x = element_text(size=12,angle=90))

f=cor.test(temp$logMut,temp$Entropy,use='pairwise.complete.obs')
f['estimate']
f['p.value']

ggplot(temp, aes(x = rscale2(logMut), y = rscale2(Entropy))) + geom_point(aes(color=binMut)) + geom_smooth(method='lm',se=F) + theme_classic() + theme(axis.line=element_line(size=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
xlab("Mutation Burden")+ylab("Transcriptional Entropy (Scaled)")+labs(color="Mutation Burden")+scale_x_continuous(expand=c(0,0),limits=c(0,1.1)) + scale_y_continuous(expand=c(0,0),limits=c(0,1.1)) + scale_color_manual(values=c('#F52F57','#A20021'))
