library(gdata)
source("./Source/toload.R")
Bdata <- list()
Bdata$id <- read.xls("portable.xls", sheet = 1, header = TRUE)
Bdata$sts <- read.xls("portable.xls", sheet = 2)
Bdata$ace <- read.xls("portable.xls", sheet = 3)
Bdata$gltB <- read.xls("portable.xls", sheet = 4)
Bdata$gmhD <- read.xls("portable.xls", sheet = 5)
Bdata$lepA <- read.xls("portable.xls", sheet = 6)
Bdata$lipA <- read.xls("portable.xls", sheet = 7)
Bdata$narK <- read.xls("portable.xls", sheet = 8)

ids <- sort(c(4, 27, 24, 845, 51, 5, 9, 20, 61, 112, 2, 6, 1, 30, 29,13, 12, 18, 8, 17, 47, 32, 35, 34, 23, 68, 21, 149, 136, 137, 138, 140, 141, 145, 63))
subBdata <- subset(Bdata$id, id %in% ids, select=c(id, st, Species))
sts <- subset(Bdata$sts, ref %in% subBdata$st)
names(sts)[1] <- "st"
subBdata <- merge(subBdata, sts, by="st", sort = FALSE)

Bseq <- rep("", length(ids))
for(j in 3:8){
    jid <- subBdata[,(j+1)]
    for(i in 1:length(ids)){
        Bseq[i] <- paste(Bseq[i],Bdata[[j]]$sequence[jid[i]], sep="")
    }
}

Bseq <- as.list(tolower(Bseq))
BList <- list()
for(i in 1:length(ids)){
    x <- Bseq[[i]]
    xDNAbin <-substring(x, seq(1,nchar(x)),seq(1,nchar(x)))
    BList[[i]] <- xDNAbin
}

speciesNames <- sapply(strsplit(as.character(subBdata$Species), split= " ", fixed=TRUE), function(x) x[2])
tipNames <- paste("B.", speciesNames, " st", subBdata$st, sep="")
names(BList) <- tipNames
names(BList) <- subBdata$st
BDNAbin <- as.DNAbin(BList)
distMat <- dist.dna(BDNAbin, as.matrix = TRUE)
initB <- EMinit(distMat, 2, diffP = 0.9)
BPhy <- as.phyDat(BDNAbin)
BResult <- DNACluster(BPhy, initP=initB)

##################################################
## sd of theta
##################################################
sdtheta <- sqrt(diag(varTheta(BResult, BPhy)))
## 95% CI of log-scale for theta:
sdlog <- sdtheta * 1/c(BResult$theta, BResult$p[1])
paste("Theta 1:","(", BResult$theta[1] * exp(-sdlog[1]),"," , BResult$theta[1] * exp(sdlog[1]),")")
paste("Theta 2:","(", BResult$theta[2] * exp(-sdlog[2]),"," , BResult$theta[2] * exp(sdlog[2]),")")
paste("P 1:","(", BResult$p[1] * exp(-sdlog[3]),"," , BResult$p[1] * exp(sdlog[3]),")")
paste("P 2:","(", BResult$p[2] * exp(-sdlog[3]),"," , BResult$p[2] * exp(sdlog[3]),")")
##################################################
## Cross validation
##################################################

initBList <- list()
for(i in 1:4){
    initBList[[i]] <- EMinit(distMat, i, diffP = 0.9)
}

cvllobs3 <- cvCluster(list(obsSeq = BPhy), 33, maxCluster = 4, mCV=100, beta=0.5, initPList = initBList)

##################################################
## Make a tree with minimun evolution
##################################################


treeUr <- fastme.bal(dist.dna(BDNAbin))
treeR <- root(treeUr,"B.oklahomensis st81")
pdf("treeB.pdf")
plot(treeR)
add.scale.bar(y = 0.5, length = 0.01)
dev.off()
## Bootstrap
boot <- boot.phylo(BDNAbin, x=dist.dna(BDNAbin), B = 2, FUN = function(xx) fastme.bal(dist.dna(xx)))
