# load HRD scores
hrdtab <- read.table("ITH_HRDscores_PC.txt",sep="\t",head=T)

# identify which patients had mixed HRD status
hrdtab$Patient <- sapply(as.character(hrdtab$Sample_Name),function(x)substr(x,start=which(strsplit(x,split="")[[1]]=="T")[1],stop=which(strsplit(x,split="")[[1]]=="T")[1]+6))
hrdtab$gt42 <- hrdtab[,'HRD_SUM']>=42

isMixed <- sapply(unique(hrdtab$Patient),function(x)length(table(hrdtab[which(hrdtab$Patient==x & hrdtab$HRD_SUM>0),'gt42'])))==2
filteredScores <- hrdtab[which(hrdtab$HRD_SUM>0),'HRD_SUM']
filteredOutcomes <- unlist(sapply(1:length(isMixed),function(x)rep(isMixed[x],length(which(hrdtab$Patient==unique(hrdtab$Patient)[x] & hrdtab$HRD_SUM>0)))))

# check the probability of mixed HRD status given HRDscore from each individual sample 
predFit <- glm(y~x,data=data.frame(y=as.numeric(filteredOutcomes),x=abs(42-filteredScores)),family="binomial")
summary(predFit)
#Coefficients:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.86002    0.26167   3.287  0.00101 ** 
#	x           -0.06780    0.01174  -5.777 7.58e-09 ***
# => v sig fit

scoreRange <- seq(from=42-max(abs(42-filteredScores)),to=42+max(abs(42-filteredScores)),length=100)
predictedProbMixed <- predict(predFit,type="response",newdata=data.frame(x=abs(42-scoreRange)))

# make Fig 4c
plot(x=scoreRange[which(scoreRange>0)],y=predictedProbMixed[which(scoreRange>0)],xlab='HRD score',ylab='Prob Mixed',pch=16)

