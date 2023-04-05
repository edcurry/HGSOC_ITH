mergedAscatTables <- list()
eventDistMats <- list()
for(i in 1:length(fullascat)){
	        cat(paste('processing data for patient',i,'of',length(fullascat),'\n'))
        mergedAscatTables[[i]] <- mergeBreakpoints(fullascat[[i]])
	        eventDistMats[[i]] <- nEventDist(mergedAscatTables[[i]])
}

meanDist <- function(x,exclude){
	        include <- setdiff(rownames(x),exclude)
        outdist <- NA
	        if(length(include)>1){
			                outdist <- mean(x[include,include][upper.tri(x[include,include])])
	        }
	        outdist
}

normalSamples <- lapply(fullascat,function(x)x$nonaberrantarrays)
fullascat.pIDs <- sapply(fullascat,function(x)substr(x$segments$sample[1],start=1,stop=7))
meanNEvents <- rep(NA,length(fullascat.pIDs))
for(i in 1:length(fullascat.pIDs)){
	        meanNEvents[i] <- meanDist(eventDistMats[[i]],exclude=normalSamples[[i]])
}

nEvents.df <- data.frame(pID=fullascat.pIDs,meanHet=meanNEvents,relapse=clin.full[fullascat.pIDs,'relapse.status'],debulking=clin.full[fullascat.pIDs,'debulk.status'],disease.status=clin.full[fullascat.pIDs,9])

# make Fig 2b
print(ggplot(nEvents.df[nEvents.df$disease.status %in% c(0,1),],aes(factor(relapse),meanHet,fill=factor(relapse))) + geom_boxplot() + geom_jitter() + facet_wrap(vars(disease.status)))

