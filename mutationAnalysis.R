# load mutation data

somatic.mut <- read.table("somatic_HighModImact_variants.csv",sep=",",head=T)
germline.mut <- read.table("germline_highimpact_variants.csv",sep=",",head=T)

# only 4 germline mutations: check these patients' tumour evolution patterns
germline.mut[,c(1,10)]
#             SAMPLE VEP.SYMBOL
#1  T14-028      BRIP1
#2 T15-058      BRCA2
#3    T17-167      BRCA2
#4  T15-022     BRCA1
germline.df <- data.frame(pID=germline.mut[,1],gene=germline.mut[,10])

# map to samples/patients, relapse-or-primary
somatic.mut$mapID <- sapply(as.character(somatic.mut[,1]),function(x)paste(strsplit(x,split="_",fixed=T)[[1]][-1],collapse="_"))
somatic.mut$pID <- sapply(as.character(somatic.mut$mapID),function(x)strsplit(x,split="_",fixed=T)[[1]][1])
somatic.mut$samplecode <- sapply(as.character(somatic.mut$mapID),function(x)rev(strsplit(x,split="_",fixed=T)[[1]])[1])

somatic.mut$mappedSample <- sapply(1:nrow(somatic.mut),function(i)ifelse(as.character(somatic.mut$pID[i]) %in% fullascat.pIDs,grep(as.character(somatic.mut$samplecode[i]),allsamples[[which(fullascat.pIDs==as.character(somatic.mut$pID[i]))]],value=T),NA))
# manually add pID for one with primary showing mutations
somatic.mut[which(somatic.mut$pID=="T16-283"),'mappedSample'] <- allsamples[[which(fullascat.pIDs=="T16-283")]]

# create dataframe with patient annotations plus mutations
mut.genes <- as.character(unique(somatic.mut$VEP.SYMBOL))
mut.samples <- setdiff(as.character(somatic.mut$mapID),NA)
mut.matrix <- array(0,dim=c(length(mut.samples),length(mut.genes)))
rownames(mut.matrix) <- mut.samples
colnames(mut.matrix) <- mut.genes
for(i in 1:nrow(somatic.mut)){
	if(!is.na(somatic.mut[i,'mapID'])){
		mut.matrix[as.character(somatic.mut$mapID[i]),as.character(somatic.mut$VEP.SYMBOL[i])] <- somatic.mut[i,'TUMOR.AF']
	}
}

# draw ggplot2 geom_tile grid with mutations per sample
# limit plot to patients with at least one non-TP53 mutation
mutatedsamples <- mut.samples[which(rowMeans(mut.matrix[,2:8])>0)]
mutatedpIDs <- unique(somatic.mut[which(somatic.mut$mapID %in% mutatedsamples),'pID'])
allmutsamples <- unique(as.character(somatic.mut[which(somatic.mut$pID %in% mutatedpIDs),'mapID']))
mutationdf <- as.data.frame(ceiling(mut.matrix[allmutsamples,]))
mutationdf$sampleID <- allmutsamples
mutationdf$pID <- sapply(allmutsamples,function(x)as.character(somatic.mut[which(somatic.mut$mapID==x),'pID'][1])) 
mutationdf$relapse <- rep('presentation',nrow(mutationdf))
mutationdf$relapse[which(mutationdf$pID %in% as.character(pairedIDs[,2]))] <- 'relapse'
mutationdf$pID[which(mutationdf$pID %in% as.character(pairedIDs[,2]))] <- sapply(mutationdf$pID[which(mutationdf$pID %in% as.character(pairedIDs[,2]))],function(x)as.character(pairedIDs[which(as.character(pairedIDs[,2])==x),1]))

mutationdf$relapse.category <- clin.full[as.character(mutationdf$pID),'relapse.category']

# make Fig 5a
library(reshape2)
mutationdf2 <- melt(mutationdf,id.var="sampleID") 
mutationdf2$pIDcode <- as.numeric(substr(gsub(as.character(mutationdf$pID),pattern="-",replace=""),start=2,stop=6))
library(ggplot2)
print(ggplot(mutationdf2,aes(reorder(sampleID,pIDcode),variable)) + geom_tile(aes(fill=value),width=0.8, height= 0.8,colour="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(0.3,"cm"),legend.key.width =unit(0.15,"cm")) +
  scale_fill_manual(values=c("#CCCCCC","#3399FF","#CC99CC","#FF007F","#FF6666","#009900","#009999","#66CC66","#FF9933","#FFCC66","#CC0000","#66FFFF","#FF33FF","#0080FF"),na.value="grey95"))

