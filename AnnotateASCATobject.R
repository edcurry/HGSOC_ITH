# get patient IDs from ASCAT output object
fullascat.pIDs <- c(sapply(ascat.out,function(x)unique(substr(x$segments$sample,start=1,stop=7))[1]),sapply(newAscat.out,function(x)unique(substr(x$segments$sample,start=1,stop=7))))

# annotate data with patient IDs
pairedIDs <- data.frame(primary=c('T14-045',"T14-093","T14-137","T15-022","T15-043","T16-046","T16-198","T16-216","T17-096","T17-167"),relapse=c("T16-175","T16-111","T17-067","T17-007","T18-098","T16-283","T19-109","T18-268","T19-038","T19-075"))
ascat.pIDs <- sapply(ascat.out,function(x)substr(x$segments$sample[1],start=1,stop=7))
pairedIDs.withCN <- pairedIDs[which(pairedIDs$primary %in% ascat.pIDs & pairedIDs$relapse %in% ascat.pIDs),]
newpairedIDs.withCN <- pairedIDs[which(pairedIDs$primary %in% fullascat.pIDs & pairedIDs$relapse %in% fullascat.pIDs),]

fullascat.primary.pIDs <- setdiff(fullascat.pIDs,as.character(newpairedIDs.withCN[,2]))

allnormalSamples <- lapply(fullascat,function(x)unique(c(x$nonaberrantarrays,grep("SN",unique(x$segments$sample),value=T),grep("BC",unique(x$segments$sample),value=T))))
allsamples <- lapply(fullascat,function(x)unique(x$segments$sample))

