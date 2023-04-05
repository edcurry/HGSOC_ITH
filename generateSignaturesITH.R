library(NMF)
library(YAPSA)
library(flexmix)

# ensure these 3 functions know where all the correct data is...
source("helper_functions.R")
source("main_functions.R")

filterSegTable <- function(cn){
    cn<-cn[!is.na(cn$segVal),]
    segTable<-c()
    for(c in unique(cn$chromosome))
    {
    snfilt<-cn[cn$chromosome==c,]
    sn.rle<-rle(snfilt[,"segVal"])
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- snfilt$start[starts[s]]
      to <- snfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(snfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    colnames(segTableRaw)<-c("chromosome","start","end","segVal")
    segTable<-rbind(segTable,segTableRaw)
    }
    segTable
}

# for each sample, first step is making the segment tables from ASCAT output:
# column headers: "chromosome", "start", "end", "segVal"
# input list must be named!
allCNsigs <- list()
allCNnames <- c()
for(i in 1:length(fullascat)){
	thisSamples <- unique(fullascat[[i]]$segments$sample)
	for(j in 1:length(thisSamples)){
		if(!thisSamples[j] %in% unlist(allnormalSamples)){
			segTable <- data.frame(chromosome=fullascat[[i]]$segments[which(fullascat[[i]]$segments$sample==thisSamples[j]),'chr'],start=fullascat[[i]]$segments[which(fullascat[[i]]$segments$sample==thisSamples[j]),'startpos'],end=fullascat[[i]]$segments[which(fullascat[[i]]$segments$sample==thisSamples[j]),'endpos'],segVal=fullascat[[i]]$segments[which(fullascat[[i]]$segments$sample==thisSamples[j]),'nMajor']+fullascat[[i]]$segments[which(fullascat[[i]]$segments$sample==thisSamples[j]),'nMinor'])
#			segTable <- filterSegTable(segTable)
			if(length(segTable)>0){
				allCNsigs[[length(allCNsigs)+1]] <- segTable
				allCNnames <- c(allCNnames,thisSamples[j])
			}
		}
	}
}
names(allCNsigs) <- allCNnames
CN_features <- extractCopynumberFeatures(allCNsigs)

sample_by_component <- generateSampleByComponentMatrix(CN_features)
allsignatures <- quantifySignatures(sample_by_component)

# save signatures object
saveRDS(allsignatures,file='ITH_allCNsignatures.rds')
