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
MultiMatch <- subset(matchedDf, Seq == "Multiple matches")$Isolate
noMatch <- subset(matchedDf, Seq == "No match")
workingDf <- subset(matchedDf, Seq != "Multiple matches" & Seq != "No match" & nchar(Seq) < 350 & nchar(Seq) > 320)
hist(unlist(lapply(workingDf$Seq,FUN=nchar)))
sink(file = "misseq.txt")
cat(c("Multiple matches:\n", MultiMatch, "\n"))
cat(c("Isolates with no matched sequences:\n", noMatch, "\n"))
sink()
connections()
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
## Distance matrix for places
##################################################
## install.packages("maps")
siteNstate <- paste(infoDf$Site, infoDf$State)
unique(siteNstate)
write.csv(sites, file = "sites.csv")


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


debug(mtCluster)
## Test mtCluster()
phyDNA <- as.DNAbin(finalPhy)
distMat <- as.matrix(dist.dna(phyDNA))
initDataP <- EMinit(distMat, 2, diffP = 0.9)
## debug(mtCluster)
## initP[[5]] <- DNADataResult$Ez[,4:8]
mtDataResult <- mtCluster(finalPhy, obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat,  initP=initPList[[2]], maxIter = 100, tol = 1e-10)
DNADataResult <- DNACluster(finalPhy, initP=initDataP)
Ez <- mtDataResult$Ez
maxEz <- apply(Ez, 1, function(x) which(x==max(x)))

hostInitP <- matrix(c(rep(0.9,5),rep(0.1, 262), rep(0.1,5), rep(0.9, 262)),267,2)
hostResult <- hostCluster(obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat, initP = randomPMat(0.1, 2, 267), model= "normal")

save.image()
cvDataListobs <- list(obsSeq=finalPhy, obsHost=hostDistMat$obsHost, distMat = hostDistMat$distMat)


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

load("result_ct.RData")


## Comparing with single trait
idHostR <- idSizeHR <- idDNA <- idHost<- idSizeDNA <- idSizeHost<- list()
for(i in 1:7){
    idDNA[[i]] <- max.col(resultD[[i]]$Ez)
    idHostR[[i]] <-  max.col(resultHR[[i]]$Ez)
    idHost[[i]] <- max.col(resultH[[i]]$Ez)
    idSizeDNA[[i]] <- table(idDNA[[i]])
    idSizeHost[[i]] <- table(idHost[[i]])
    idSizeHR[[i]] <- table(idHost[[i]])
}

## Use normal likelihood for distance between hosts.
load("resultN.RData")

idN <- idHN <- idRN <- idHRN<- idSizeN <- idSizeHN <- idSizeRN <- idSizeHRN<- list()
## HN and HRN do not cluster?
for(i in 1:7){
    idN[[i]] <- max.col(resultN[[i]]$Ez)
    idHN[[i]] <-  max.col(resultHN[[i]]$Ez)
    idRN[[i]] <- max.col(resultRN[[i]]$Ez)
    idHRN[[i]] <- max.col(resultHRN[[i]]$Ez)
    idSizeN[[i]] <- table(idN[[i]])
    idSizeHN[[i]] <- table(idHN[[i]])
    idSizeRN[[i]] <- table(idRN[[i]])
    idSizeHRN[[i]] <- table(idHRN[[i]])
}

apply(cvllobsN, 2, mean)
apply(cvllobsRN, 2, mean)

## table choosing the number of clusters:
cvTable <- matrix(, 3, 7)
colnames(cvTable) <- as.character(1:7)
rownames(cvTable) <- c("log likelihood","cv log likelihood","BIC")
for( i in 1:7){
    cvTable[1,i] <- resultRN[[i]]$logLike
    cvTable[2,i] <- mean(cvllobsN[,i])
    cvTable[3,i] <- 2 * (log(267) * (267 + 34 +2) * i - resultRN[[i]]$logLike)
}
length(cvllobsCTN)
xtable(cvTable, digits=0)

idSizeDNA[[4]]
idSize[[4]]

## Host distribution:
host <- list()
for(i in 1:4){
    host[[i]] <- hostDistMat$obsHost[which(identification[[4]]==i)]
}
## 18 6 isolates in 2,3,4
## 23 in 1(3), 2, 3(mostly), 4
## 22 in 4(5),1(15) 1(1)
## 17 in 2(mostly), 1(3)
## 25 in 1(32), 2(3)
## 35
## 20 in 2 adn 3

## Changes
idSize[[4]]
idSizeDNA[[4]]
hostDNA <- list()
for(i in 1:4){
    hostDNA[[i]] <- hostDistMat$obsHost[which(idDNA[[4]]==i)]
}

sort(hostDNA[[4]])
sort(host[[2]])

sort(hostDNA[[2]])
sort(host[[1]])
##  29 split into 29(3) and 29(15)


## pie chart:
slices <- lbls <- colorPie <- list()
hostNames <- paste(hostDf[,2], hostDf[,1])
for(i in 1:2){
    slices[[i]] <- 0
    lbls[[i]] <- 0
    unisp <- unique(host[[i]])
    colorPie[[i]] <- unisp
    for(j in 1:length(unisp)){
        slices[[i]][j] <- sum(host[[i]] == unisp[j])
        lbls[[i]][j] <- paste(hostNames[unisp[j]], "\n", slices[[i]][j], sep="")
    }
}

pie(slices[[1]], lbls[[1]])
pie(slices[[2]], lbls[[2]])

a <- c(2:4,6:10)
b <- unique(c(colorPie[[3]],colorPie[[4]]))
for(i in 1:length(b)){
    colorPie[[3]][which(colorPie[[3]] == b[i])] <- a[i]
    colorPie[[4]][which(colorPie[[4]] == b[i])] <- a[i]
}

pie(slices[[3]], lbls[[3]], col=colorPie[[3]])
pie(slices[[4]], lbls[[4]],col=colorPie[[4]])



 mtResult <- mtCluster(finalPhy, obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat,  initP=initPList[[2]], maxIter = 100, tol = 1e-6)
DNAresult <-  DNACluster(finalPhy,  initP=initPList[[2]], maxIter = 100, tol = 1e-6)

idDNA <- max.col(mtResult$Ez)
idmt <-  max.col(DNAresult$Ez)

## Host distribution:
hostDNA <- hostMt<- list()
for(i in 1:2){
    hostDNA[[i]] <- sort(hostDistMat$obsHost[which(idDNA==i)])
    hostMt[[i]] <- sort(hostDistMat$obsHost[which(idmt==i)])
}
## pie chart:
## DNA
slices <- lbls <- colorPie <- list()
hostNames <- paste(hostDf[,2], hostDf[,1])
for(i in 1:2){
    slices[[i]] <- 0
    lbls[[i]] <- 0
    unisp <- unique(hostDNA[[i]])
    colorPie[[i]] <- unisp
    for(j in 1:length(unisp)){
        slices[[i]][j] <- sum(hostDNA[[i]] == unisp[j])
        lbls[[i]][j] <- paste(hostNames[unisp[j]], "\n", slices[[i]][j], sep="")
    }
}

pie(slices[[1]], lbls[[1]])
pie(slices[[2]], lbls[[2]])

## pie chart:
## mt
slices <- lbls <- colorPie <- list()
hostNames <- hostDf[,1]
for(i in 1:2){
    slices[[i]] <- 0
    lbls[[i]] <- 0
    unisp <- unique(hostMt[[i]])
    colorPie[[i]] <- unisp
    for(j in 1:length(unisp)){
        slices[[i]][j] <- sum(hostMt[[i]] == unisp[j])
        lbls[[i]][j] <- paste(hostNames[unisp[j]], "\n", slices[[i]][j], sep="")
    }
}

pie(slices[[1]], lbls[[1]])
pie(slices[[2]], lbls[[2]])
unisp
