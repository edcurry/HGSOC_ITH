clin.full <- read.table("ITHclinicaldata.csv",sep="\t",head=T)
rownames(clin.full) <- as.character(clin.full[,2])
clin.full$Relapsed <- as.numeric(sapply(as.character(clin.full$Date.of.relapse),function(x)!x %in% c("",NA,"?")))
clin.full$Dead <- as.numeric(sapply(as.character(clin.full$Date.of.death),function(x)!x %in% c("",NA,"?")))
clin.full$PFS.months <- sapply(as.character(clin.full[,"PFS.from.surgery.date"]),function(x)as.numeric(strsplit(x,split=" ")[[1]][1]))
clin.full$relapse.status <- rep(NA,nrow(clin.full))
clin.full$relapse.status[which(clin.full$Relapsed==0 & clin.full$PFS.months>12)] <- 'no.relapse'
clin.full$relapse.status[which(clin.full$Relapsed==1 & clin.full$PFS.months>=12)] <- 'relapse>1yr'
clin.full$relapse.status[which(clin.full$Relapsed==1 & clin.full$PFS.months<12)] <- 'relapse<1yr'
clin.full$relapse.status[which(clin.full$Relapsed==1 & clin.full$PFS.months<6)] <- 'relapse<6mo'

clin.full$debulk.status <- c('residual','none')[1+as.numeric(toupper(clin.full[,'Optimal.debulking.'])=="YES")]

clin.full$relapse.category <- NA
clin.full[which(clin.full$relapse.status=='no.relapse'),'relapse.category'] <- 'x.no.relapse'
clin.full[which(clin.full$relapse.status=='relapse>1yr'),'relapse.category'] <- 'sensitive'
clin.full[which(clin.full$relapse.status %in% c('relapse<1yr','relapse<6mo')),'relapse.category'] <- 'resistant'
pIDs <- unique(rownames(clin.full))

# get median time-to-last-follow-up
clin.full$Follow.up.time <- as.Date(clin.full$Date.of.last.follow.up,format="%d/%m/%Y")-as.Date(clin.full$Date.of.surgery,format="%d/%m/%Y")
12*(as.numeric(median(clin.full$Follow.up.time))/365)

