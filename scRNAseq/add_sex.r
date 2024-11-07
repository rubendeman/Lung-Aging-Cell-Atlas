library(readxl)
BCM_controls_metadata <- read_excel("BCM_controls_metadata.xlsx")
BCM_controls_metadata <- BCM_controls_metadata[,c(3,6,7)]

ipf <- read.csv('sc_demo.csv')

BCM_controls_metadata <- rbind(BCM_controls_metadata,setNames(ipf,names(BCM_controls_metadata)))
#BCM_controls_metadata$Gender[BCM_controls_metadata$Gender=='Female']<-0
#BCM_controls_metadata$Gender[BCM_controls_metadata$Gender=='Male']<-1

sex <- BCM_controls_metadata$Gender[match(immune.combined@meta.data$orig.ident,as.character(unlist(BCM_controls_metadata[,1])))]
immune.combined<-AddMetaData(immune.combined,metadata=as.factor(sex),col.name='sex')