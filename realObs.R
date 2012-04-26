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
library(phangorn)

finalList <- list()
for(i in 1:length(finalDf$Seq)){
    x <- finalDf$Seq[i]
    xDNAbin <-substring(x, seq(1,nchar(x)),seq(1,nchar(x)))
    finalList[[i]] <- xDNAbin
}

names(finalList) <- finalDf$Isolate
finalDNAbin <- as.DNAbin(finalList, fill.with.gaps = TRUE)


## Align sequences data
finalAlign <- clustal(finalDNAbin)
image.DNAbin(finalAlign)
finalPhy <- as.phyDat(finalAlign)
dm <- dist.dna(finalDNAbin)


## save
save(finalPhy, file = "align_data.RData")
unorderPhy <- finalPhy
for(i in 1:length(unorderPhy)){
    j <- which(names(unorderPhy)[i] == names(finalList))
    finalPhy[[j]] <- unorderPhy[[i]]
    names(finalPhy)[j] <- names(unorderPhy)[i]
}

##################################################
## Clustering
##################################################

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


## run script.R
## after running script.R
load("result.RData")

apply(cvllobs, 2, sd)
apply(cvllobs, 2, mean)

result[[20]]$p

identification <- idSize <- list()
for(i in 1:7){
    identification[[i]] <- max.col(result[[i]]$Ez)
    idSize[[i]] <- table(identification[[i]])
}

idSize[[7]]

g1 <- which(identification[[20]] == 20) # 62
g4 <- which(identification[[20]] == 15) # 56
g6 <- which(identification[[20]] == 3) # 35 20[3] + 20[15] == 9[3]
g2 <- which(identification[[20]] == 18) # 31 == to 15[8] + 15[13] == 3[1]
g3 <- which(identification[[15]] == 4) # 14 == to 20[14]+[13]
g5 <- which(identification[[20]] == 10) # 38
g7 <-  which(identification[[15]] == 3) # unassigned: [[20]][6](11) [12](9) [17](11) == [[15]][3](31)

## g4+g6 <- 8[7] 91
## g5+g7 == 8[8] 69

nrP <- function(group,n){
    M <- length(group)
    PMat <- matrix(, n, M)
    TC <- 1.1/ M
    WC <- (1 - TC) / M
    for(i in 1:M){
        PMat[group[[i]],i] <- TC
        PMat[-group[[i]],] <- WC
    }
    PMat
}

glist5 <- list(g1,g2,g3,c(g4,g6), c(g5,g7))
glist6 <- list(g1,g2,g3,g4,g6, c(g5,g7))
glist7 <- list(g1,g2,g3,g4,g6, g5,g7)

set.seed(42)
initP <- list()
initP[[1]] <- matrix(1, 267, 1)
initP[[2]] <- randomPMat(0.01, 2, 267)
initP[[3]]<- randomPMat(0.01, 3, 267)
initP[[4]]<- randomPMat(0.01, 4, 267)
initP[[5]] <- nrP(glist5,267)
initP[[6]] <- nrP(glist6,267)
initP[[7]] <- nrP(glist7,267)
