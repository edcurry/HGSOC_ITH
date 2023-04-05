library(NMF)

# read in pre-computed CN signature profiles
allsignatures <- readRDS('ITH_allCNsignatures.rds')

# for each patient (primary samples only), compute euclidean distance of sig scores
CNShet <- function(x){
        dmat <- as.matrix(dist(t(x)))
        mean(dmat[upper.tri(dmat)],na.rm=T)
} 

allprimary.sigmats <- list()
for(i in 1:length(fullascat.primary.pIDs)){
	this.samp <- grep(fullascat.primary.pIDs[i],colnames(allsignatures),value=T)
	allprimary.sigmats[[i]] <- allsignatures[,this.samp]
}
allprimary.sighet <- sapply(allprimary.sigmats,CNShet)
nsamps <- sapply(allprimary.sigmats,ncol)

# make Fig 3b, requires clinical data
require(ggplot2)
CNsig.df <- data.frame(pID=fullascat.primary.pIDs[which(nsamps==5)],CNsig.het=allprimary.sighet[which(nsamps==5)],relapse=clin.full[fullascat.primary.pIDs[which(nsamps==5)],'relapse.category'])
print(ggplot(CNsig.df,aes(x=factor(relapse),y=CNsig.het)) + geom_boxplot(aes(fill=factor(relapse))) + geom_jitter(aes(color=factor(relapse))))

# now compare NEvent distances to CNsig distances (requires sequential running of steps in README)
# make Fig 3d
allprimary.eventhet <- sapply(fulleventdist.bypatient,function(x)mean(x[upper.tri(x)],na.rm=T))
plot(x=allprimary.eventhet[which(nsamps.event==5)],y=allprimary.sighet[which(nsamps==5)],xlab='NEvent Dist',ylab='CNsig Dist')

# draw example signature profiles for a patient, as in Fig 3a
require(ggplot2)
t18085.sigdf <- data.frame(score=as.numeric(allprimary.sigmats[[which(fullascat.primary.pIDs=='T18-085')]]),CNsig=rep(1:7,ncol(allprimary.sigmats[[which(fullascat.primary.pIDs=='T18-085')]])),sampleNo=rep(colnames(allprimary.sigmats[[which(fullascat.primary.pIDs=='T18-085')]]),each=7))
print(ggplot(t18085.sigdf,aes(factor(sampleNo),score,fill=factor(CNsig))) + geom_col())

t15085.sigdf <- data.frame(score=as.numeric(allprimary.sigmats[[which(fullascat.primary.pIDs=='T15-085')]]),CNsig=rep(1:7,ncol(allprimary.sigmats[[which(fullascat.primary.pIDs=='T15-085')]])),sampleNo=rep(colnames(allprimary.sigmats[[which(fullascat.primary.pIDs=='T15-085')]]),each=7))
print(ggplot(t15085.sigdf,aes(factor(sampleNo),score,fill=factor(CNsig))) + geom_col())

# compute patient-wise prim-relapse matrices
CNShet <- function(x){
	dmat <- as.matrix(dist(t(x)))
        mean(dmat[upper.tri(dmat)],na.rm=T)
} 

allpaired.sigmats <- list()
for(i in 1:nrow(pairedIDs)){
	allpaired.sigmats[[i]] <- list()
	this.primsamp <- grep(as.character(pairedIDs[i,1]),colnames(allsignatures),value=T)
        allpaired.sigmats[[i]]$primary <- allsignatures[,this.primsamp]
	this.relsamp <- grep(as.character(pairedIDs[i,2]),colnames(allsignatures),value=T)
	allpaired.sigmats[[i]]$relapse <- allsignatures[,this.relsamp]
}

allpaired.primsighet <- sapply(allpaired.sigmats,function(x)CNShet(x$primary))

# check method for computing distances:
sqrt(sum((allpaired.sigmats[[1]]$primary[,1]-allpaired.sigmats[[1]]$primary[,2])^2))
#[1] 0.1202153
as.matrix(dist(t(allpaired.sigmats[[1]]$primary)))[1,2]
#[1] 0.1202153

allpaired.sighet <- rep(NA,length(allpaired.sigmats))
for(i in 1:length(allpaired.sighet)){
	if(is.null(ncol(allpaired.sigmats[[i]]$relapse))){
		allpaired.sighet[i] <- mean(apply(allpaired.sigmats[[1]]$primary-allpaired.sigmats[[1]]$relapse,2,function(x)sqrt(sum(x^2))),na.rm=T) 
	}
	if(!is.null(ncol(allpaired.sigmats[[i]]$relapse))){
		if(ncol(allpaired.sigmats[[i]]$relapse)>0){
			thisDmat <- array(NA,dim=c(ncol(allpaired.sigmats[[i]]$primary),ncol(allpaired.sigmats[[i]]$relapse)))
			for(j in 1:nrow(thisDmat)){
				for(k in 1:ncol(thisDmat)){
					thisDmat[j,k] <- sqrt(sum((allpaired.sigmats[[i]]$primary[,1]-allpaired.sigmats[[i]]$relapse[,1])^2))
				}
			}
			allpaired.sighet[i] <- mean(thisDmat,na.rm=T)
		}
	}
}

