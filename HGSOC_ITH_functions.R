# for CNA profiles, generate merged-breakpoints binary profile for each sample

require(GenomicRanges)

mergeBreakpoints <- function(x){
	allchrs <- unique(x$segments$chr)
	allbkps <- lapply(allchrs,function(y)unique(sort(c(x$segments$startpos[x$segments$chr==y],x$segments$endpos[x$segments$chr==y]))))
	allsamples <- unique(x$segments$sample)
	allbkpt.grs <- list()
	for(i in 1:length(allchrs)){
		allbkpt.grs[[i]] <- GRanges(seqnames=allchrs[i],ranges=IRanges(start=allbkps[[i]][1:(length(allbkps[[i]])-1)],end=allbkps[[i]][2:length(allbkps[[i]])]-1))
	}
	allsample.grs <- lapply(allsamples,function(y)GRanges(seqnames=x$segments$chr[x$segments$sample==y],ranges=IRanges(start=x$segments$startpos[x$segments$sample==y],end=x$segments$endpos[x$segments$sample==y]),nMajor=x$segments$nMajor[x$segments$sample==y],nMinor=x$segments$nMinor[x$segments$sample==y]))
	merged.grs <- list()
	for(i in 1:length(allsamples)){
		thismerged.grs <- list()
		for(j in 1:length(allbkpt.grs)){
			thismerged.grs[[j]] <- allbkpt.grs[[j]]
			#thismerged.grs[[j]]$nMajor <- sapply(allbkpt.grs[[j]],function(y)mean(subsetByOverlaps(allsample.grs[[i]],y)$nMajor))
			#thismerged.grs[[j]]$nMinor <- sapply(allbkpt.grs[[j]],function(y)mean(subsetByOverlaps(allsample.grs[[i]],y)$nMinor))
			# should prob assume nMajor=1 and nMinor=1 for all implicit regions
			thismerged.grs[[j]]$nMajor <- rep(NA,length(allbkpt.grs[[j]]))
			thismerged.grs[[j]]$nMinor <- rep(NA,length(allbkpt.grs[[j]]))
			for(k in 1:length(allbkpt.grs[[j]])){
				thismerged.grs[[j]]$nMajor[k] <- mean(subsetByOverlaps(allsample.grs[[i]],allbkpt.grs[[j]][k])$nMajor)
				thismerged.grs[[j]]$nMinor[k] <- mean(subsetByOverlaps(allsample.grs[[i]],allbkpt.grs[[j]][k])$nMinor)
			}
		}
		merged.grs[[i]] <- thismerged.grs[[1]]
		for(j in 1:length(allbkpt.grs)){
			merged.grs[[i]] <- c(merged.grs[[i]],thismerged.grs[[j]])
		}
	}
	merged.df <- data.frame(sample=allsamples[1],chr=seqnames(merged.grs[[1]]),startpos=start(merged.grs[[1]]),endpos=end(merged.grs[[1]]),nMajor=merged.grs[[1]]$nMajor,nMinor=merged.grs[[1]]$nMinor)
	for(i in 2:length(allsamples)){
		merged.df <- rbind(merged.df,data.frame(sample=allsamples[i],chr=seqnames(merged.grs[[i]]),startpos=start(merged.grs[[i]]),endpos=end(merged.grs[[i]]),nMajor=merged.grs[[i]]$nMajor,nMinor=merged.grs[[i]]$nMinor))
	}
	merged.df
}

# use this to find number of discriminatory events between each pair

nEventDist <- function(x){
	allsamples <- unique(x$sample)
	allchrs <- unique(x$chr)
	Dmat <- array(0,dim=rep(length(allsamples),2))
	rownames(Dmat) <- allsamples
	colnames(Dmat) <- allsamples
	for(i in 1:length(allsamples)){
		profile1 <- x[which(x$sample==allsamples[i]),]
		for(j in 1:length(allsamples)){
			if(j!=i){
				profile2 <- x[which(x$sample==allsamples[j]),]
				Dmat[i,j] <- length(unique(c(which((profile1$nMajor-profile2$nMajor)!=0),which((profile1$nMinor-profile2$nMinor)!=0))))
			}
		}
	}
	Dmat
}

getTumourDists <- function(x,exclude){
        xt <- x[setdiff(rownames(x),exclude),setdiff(colnames(x),exclude)]
        xt[upper.tri(xt)]
}

classifySamples <- function(x,exclude,cutoff){
	        xt <- x[setdiff(rownames(x),exclude),setdiff(colnames(x),exclude)]
        xp <- c(grep('vary',colnames(xt),value=T),grep('rimary',colnames(xt),value=T))
	        class <- rep(NA,ncol(xt))
	        if(length(xp)==1) class <- as.numeric(xt[xp,]>cutoff)
		        if(length(xp)>1){
				                class <- as.numeric(colMeans(xt[xp,])>cutoff)
		        }
		        class
}

getClones <- function(x,exclude,cutoff){
	        xt <- x[setdiff(rownames(x),exclude),setdiff(colnames(x),exclude)]
        clones <- 1
	        if(length(xt)>1) clones <- cutree(hclust(as.dist(xt)),h=cutoff)
	        clones
}

colourByClone <- function(x){
	  if (is.leaf(x)) {
		      ## fetch label
		      label <- attr(x, "label")
#cat(paste("label:",label,"\n"))
    labelClass <- as.numeric(substr(label,start=1,stop=1))
        attr(x, "edgePar") <- c(attributes(x)$edgePar,list(col=c('blue','green','red','yellow','purple','magenta')[labelClass]))
      }
  x
}
