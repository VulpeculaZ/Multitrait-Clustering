###########################################################
## Process read sequence and host data
##################################################


###########################################################
## Read sequence data and clean up
##################################################
source("./Source/readSeq.R")
source("./Source/DNAcluster.R")
## read sequence file from a directory to a data frame.
seqDf <- readSeq(path="./seqs_new")
infoDf <- read.csv("seqinfo.csv", stringsAsFactors =FALSE, sep=";")
matchedDf <- matchSeq(seqDf, infoDf)
workingDf <- subset(matchedDf, Seq != "Multiple matches" & Seq != "No match" & nchar(Seq) < 350 & nchar(Seq) > 320)
hist(unlist(lapply(workingDf$Seq,FUN=nchar)))


##################################################
## Distance matrix for host species
##################################################
library(openNLP)

hostDf <- read.csv(file = "rodent taxonomy.csv", sep =";")
## substr(hostDf$Species, regexpr(' ', hostDf$Species)+1, stop=nchar(as.character(hostDf$Species)))
hostDf$Species <- sapply(X = hostDf$Species, FUN = tokenize)[2,]
hostDf <- subset(hostDf,select = c(Species, Genus, Subfamily., Family, Suborder))
## Change the name to only species name
workingDf$Species <- substr(workingDf$Species, regexpr('[\\.[:space:]]+', workingDf$Species) + 1, stop=nchar(workingDf$Species))

## Clean up and correct typo
reg <- 'ordi$'
workingDf$Species <- sub(reg, 'ordii', workingDf$Species)
workingDf$Species <- sub('boylei', 'boylii', workingDf$Species)
excluded <- c('nuttalli', 'humulis', 'megalotis', 'bottae', 'quadrivittatus', 'musculus')
finalDf <- workingDf[!is.element(workingDf$Species, excluded),]

## Calculate distance matrix
hostDistMat <- hostDist(hostDf, finalDf)
save(hostDistMat, file = "host.RData")

##################################################
## Alignment
##################################################
## Prepare for alignment
finalMat <- substring(finalDf$Seq[[1]], seq(1,nchar(finalDf$Seq[[1]])), seq(1,nchar(finalDf$Seq[[1]])))
for(i in 2:length(finalDf$Seq)){
    x <- finalDf$Seq[i]
    xDNAbin <-substring(x, seq(1,nchar(x)),seq(1,nchar(finalDf$Seq[[i]])))
    finalMat <- rbind(finalMat, xDNAbin)
}
rownames(finalMat) <- finalDf$Isolate
finalDNAbin <- as.DNAbin(finalMat, fill.with.gaps = TRUE)


finalList <- list()
for(i in 1:length(finalDf$Seq)){
    finalList[[i]] <- substring(finalDf$Seq[[i]], seq(1,nchar(finalDf$Seq[[i]])), seq(1,nchar(finalDf$Seq[[i]])))
}
names(finalList) <- finalDf$Isolate
finalDNAbin <- as.DNAbin(finalList, fill.with.gaps = TRUE)

## Align sequences data
finalAlign <- clustal(finalDNAbin)
image.DNAbin(finalAlign)
finalPhy <- as.phyDat(finalAlign)

## save
save(finalPhy, hostDistMat, file = "obs_data.RData")

##################################################
## Clustering
##################################################

library(phangorn)
## library(abind)
library(inline)
library(RcppArmadillo)
## library(nnls)

load("obs_data.RData")
load("host.RData")
postRcpp <- cxxfunction(signature(data ="list",  P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/likelihood.cpp"), collapse = "\n" ))
cvRcpp <- cxxfunction(signature(data ="list", cvData="list",P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/cvlikelihood.cpp"), collapse = "\n" ))

randomPMat <- function(percentDiff, m, n){
    p <- runif(m*n, (1-percentDiff)/n, (1+percentDiff)/n)
    pMat <- matrix(p, n, m)
    pMat <- pMat / rowSums(pMat)
}

initDataP <- randomPMat(0.1, 6, 267)

debug(mtCluster)

mtDataResult <- mtCluster(finalPhy, P = P, obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat, M = 6, initP=initDataP, maxIter = 100, tol = 1e-10, CV=list(obsHost=hostDistMat$obsHost, seq=finalPhy))
Ez <- mtDataResult$Ez
maxEz <- apply(Ez, 1, function(x) which(x==max(x)))

save.image()
cvDataListobs <- list(obsSeq=finalPhy, P=P, obsHost=hostDistMat$obsHost, distMat = hostDistMat$distMat)

## 2: 8:15

## 3: 15:20
cvllobs3 <- cvCluster(cvDataListobs, 267, maxCluster = c(15:20), mCV=100, beta=0.5)
save.image()


cvSimObs <- rep(0,3)
for(i in 1:100){
    cvSimObs <- cvSimObs + (cvllobs[[i]] == max(cvllobs[[i]]))
}

## the lowest cv likelihood is at 15th clusters
