fullascat.primary <- fullascat[which(fullascat.pIDs %in% fullascat.primary.pIDs)]
interpatient.dists <- rep(NA,500)
for(i in 1:length(interpatient.dists)){
	        if(i %% 10 == 0) cat(paste('processing pair',i,'\n'))
        p <- sample(length(fulleventdist.bypatient),2,replace=F)
	        s1 <- sample(setdiff(colnames(fulleventdist.bypatient[[p[1]]]),unlist(allnormalSamples)),1)
	        s2 <- sample(setdiff(colnames(fulleventdist.bypatient[[p[2]]]),unlist(allnormalSamples)),1)
		        mergedAscat <- list()
		        mergedAscat$segments <- fullascat.primary[[p[1]]]$segments[which(fullascat.primary[[p[1]]]$segments$sample==s1),]
			        mergedAscat$segments <- rbind(mergedAscat$segments,fullascat.primary[[p[2]]]$segments[which(fullascat.primary[[p[2]]]$segments$sample==s2),])
			        ttab <- mergeBreakpoints(mergedAscat)
				        tmat <- nEventDist(ttab)
				        interpatient.dists[i] <- mean(getTumourDists(tmat,exclude=unlist(allnormalSamples)))
}

intrapatient.dists <- list()
for(i in 1:length(fulleventdist.bypatient)){
	        cat(paste('processing patient',i,'\n'))
        this.samples <- setdiff(colnames(fulleventdist.bypatient[[i]]),unlist(allnormalSamples))
	        intrapatient.dists[[i]] <- NA
	        for(j in 1:(length(this.samples)-1)){
			                for(k in (j+1):length(this.samples)){
						                        mergedAscat <- fullascat.primary[[i]]
		                        mergedAscat$segments <- mergedAscat$segments[which(mergedAscat$segments$sample %in% this.samples[c(j,k)]),]
					                        ttab <- mergeBreakpoints(mergedAscat)
					                        tmat <- nEventDist(ttab)
								                        intrapatient.dists[[i]] <- c(intrapatient.dists[[i]],mean(getTumourDists(tmat,exclude=unlist(allnormalSamples))))
								                }
		        }
		        intrapatient.dists[[i]] <- setdiff(intrapatient.dists[[i]],NA)
}
intrapatient.dists.individual <- unlist(intrapatient.dists)

dist.table <- data.frame(x=c(intrapatient.dists.individual,interpatient.dists),y=c(rep(0,length(intrapatient.dists.individual)),rep(1,length(interpatient.dists))))
dist.fit <- glm(y~x,data=dist.table,family='binomial')
# find cutoff where prediction flips to inter-patient
dist.range <- seq(from=1,to=max(dist.table$x),length=10000)
dist.predout <- predict(dist.fit,newdata=data.frame(x=dist.range),type='response')
dist.cutoff <- dist.range[which(dist.predout>0.5)[1]]

# make Fig 2d
require(ggplot2)
dist.plotdf <- data.frame(nEvents=dist.table$x,patient=c('Same-patient','Other')[1+dist.table$y])
print(ggplot(dist.plotdf,aes(nEvents,fill=patient,colour=patient)) + geom_density(alpha=0.4))

nClones <- sapply(fulleventdist.bypatient[which(nsamps.event==5)],function(x)length(unique(getClones(x,exclude=unlist(allnormalSamples),cutoff=dist.cutoff))))

nClones.df <- data.frame(pID=names(fulleventdist.bypatient)[which(nsamps.event==5)],nClones=nClones,relapse=clin.full[names(fulleventdist.bypatient)[which(nsamps.event==5)],'relapse.category'])

# make Fig 2e
print(ggplot(nClones.df,aes(relapse,nClones,fill=relapse)) + geom_boxplot() + geom_jitter())

# make clone-labelled dendrograms for primary samples only as in Fig 2f
fullpatients.nclones <- rep(NA,length(fullascat.primary.pIDs))
for(i in 1:length(fullascat.primary.pIDs)){
        if(length(which(substr(rownames(fulleventdists),start=1,stop=7)==fullascat.primary.pIDs[i]))>1){
                this.mat <- fulleventdists[substr(rownames(fulleventdists),start=1,stop=7)==fullascat.primary.pIDs[i],substr(colnames(fulleventdists),start=1,stop=7)==fullascat.primary.pIDs[i]]
                this.clones <- getClones(this.mat,exclude=c(unlist(normalSamples),unlist(newNormals)),cutoff=dist.cutoff)
                fullpatients.nclones[i] <- length(unique(this.clones))
                colnames(this.mat) <- paste(this.clones,colnames(this.mat),sep=":")
                rownames(this.mat) <- paste(this.clones,rownames(this.mat),sep=":")
                this.dendro <- as.dendrogram(hclust(as.dist(this.mat)))
                coloured.dendro <- dendrapply(this.dendro,colourByClone)
                png(file=paste('ITH_',fullascat.primary.pIDs[i],'_EventDistHclust_primary_labelled.png',sep=''),width=6,height=6,units='in',res=300)
                par(mar=c(10,4,4,2)+0.1)
                plot(coloured.dendro,ylab='nCNA')
                dev.off()
        }
}

# make clone-labelled dendrograms for primary and relapse samples as in Fig 2g
for(i in 1:length(pairedDistMats)){
        this.mat <- pairedDistMats[[i]][!rownames(pairedDistMats[[i]]) %in% unlist(normalSamples),!colnames(pairedDistMats[[i]]) %in% unlist(normalSamples)]
        this.clones <- getClones(this.mat,exclude=unlist(normalSamples),cutoff=dist.cutoff)
        colnames(this.mat) <- paste(this.clones,colnames(this.mat),sep=":")
        rownames(this.mat) <- paste(this.clones,rownames(this.mat),sep=":")
        this.dendro <- as.dendrogram(hclust(as.dist(this.mat)))
        coloured.dendro <- dendrapply(this.dendro,colourByClone)
        png(file=paste('ITH_',pairedIDs.withCN[i,'primary'],'_EventDistHclust_withRelapse_labelled.png',sep=''),width=6,height=6,units='in',res=300)
        par(mar=c(10,4,4,2)+0.1)
        plot(coloured.dendro,ylab='nCNA')
        dev.off()
}

