


##################################################
## Distance matrix for sequences
## load the libraries
source("http://www.bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)
library(ape)
data(woodmouse)

## preparation
nSeq <- dim(workingDf)[1]
distMat <- matrix(nrow=nSeq, ncol=nSeq)
rownames(distMat) <- workingDf$Isolate
colnames(distMat) <- workingDf$Isolate

## specify the substitution matrix
subMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)

masim <- sapply(1:(nSeq-1), function(x) sapply((x+1):nSeq, function(y) try(pid(pairwiseAlignment(workingDf$Seq[x], workingDf$Seq[y], substitutionMatrix=subMat, gapOpening = -5, gapExtension = -2)))))

## get the symetric distance matrix and diagonal element.
for(i in 1:(nSeq-1)){
    distMat[(i+1):nSeq, i] <- distMat[i, (i+1):nSeq] <- 1 - masim[[i]]/100
}
diag(distMat) <- 0






##################################################
## dist mat for clustering
##################################################

set.seed(42)
tb <- read.tree(text = "((a,b,c),(d,e,f));")
tb <- reorderPruning(tb)
tb$edge.length <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2)
plot(tb)

library(nnls)
## testing least square fitting on DNA sequences
sb <- simSeq(tb, l=500)
X <- designTree(tb)
X <- designUnrooted(tb)
sDNAb <- as.DNAbin(sb)
distb <- as.matrix(dist.dna(sDNAb))
y <- distb[lower.tri(distb)]
a <- tree.ls(y, X)
b <- nnls(X, y)
a <- nnls.tree(distb, tb)
b <- nnls(dataTree$X, dataTree$y)
b <- cluster.nnls(dataTree$y, dataTree$X)
index <- c(1,5)
library(Matrix)
library(MatrixModels)

dataTree <- designObsClade(distb)
initP <- matrix(c(rep(0.55,3),rep(0.45,3), rep(0.45,3), rep(0.55,3)), 6,2)

w <- combp(dataTree$X, initP)
b <- cluster.nnls(dataTree$y, dataTree$X, w = w)


set.seed(42)
b.result <- list()
for(i in 1:1000){
    p <- 0.5 + runif(12, 0, 0.1)
    initP <- matrix(p, 6, 2)
    b <- cluster.em(dataTree$y, dataTree$X, pMat = initP, maxIter = 100)
    b.result[[i]] <- b$Ez
}
nCorrect(tc = tb, r = b.result)


##################################################
## Use real distance matrix from mammalian hosts:
##################################################
mamData <- designObsClade(hostDistMat)
initMamP <- randomPMat(0.1, 6, 267)
mamR1 <- cluster.em(mamData$y, mamData$X, pMat = initMamP, maxIter = 10)

##################################################
## EM aglorithm for host with mixed MVnormal
##################################################

EMmix <- function(x, theta, fixsig = TRUE) {
    mu <- theta$mu
    sigma <- theta$sigma
    p <- theta$p
    M <- dim(mu)[1]
    uniqueX <- diag(dim(x)[2])
    lookup <- apply(x,1,function(x) which(x==1))
    ## E step
    Ez <- matrix(, dim(x)[1], M)
    dXp <- matrix(, dim(x)[2], M)
    for(j in 1:M){
        distX <- 0
        for(k in 1:dim(x)[2]){
            distX[k] <- sqrt(sum((uniqueX[k,] - mu[j,])^2))
        }
        dXp[,j] <- p[j] * dnorm(x=distX, sd=sigma[j])
        Ez[,j] <- dXp[lookup, j]
    }
    Ez <- sweep(Ez, 1, rowSums(Ez), "/")
    colSums.Ez <- colSums(Ez)
    ## M step
    mu.new <- mu
    for(j in 1:M){
        mu.new[j,] <- colSums(x * Ez[,j])
    }
    mu.new <- mu.new / colSums.Ez
    sqRes <- matrix(0, dim(x)[1], M)
    for(j in 1:M){
        sqRes[,j] <- apply(x, 1, function(x) sum((x-mu.new[j,])^2))
    }
    sigma.new <- sqrt(colSums(Ez * sqRes) / colSums.Ez)
    p.new <- colSums.Ez / sum(colSums.Ez)
    ## pack up result
    list(mu = mu.new, sigma = sigma.new, p = p.new)
}


EMmixFit <- function(x, theta, maxIter, tol){
    for( i in 1:maxIter){
        thetaNew <- EMmix(x, theta)
        if(sum( (thetaNew$mu - theta$mu)^2 ) < tol) break
        else theta <- thetaNew
    }
    theta$nIter <- i
    theta
}


testMV <- matrix(0,10,3)
testMV[1:5,1] <- 1
testMV[6:8,2] <- 1
testMV[9:10,3] <- 1
init <- list(p=c(0.5,0.5), sigma=c(0.2,0.2), mu = matrix(c(0.4,0.2,0.3,0.4,0.3,0.4),2,3))


EMmix(testMV, init, fixsig = TRUE)
EMmixFit(testMV, init, tol=1E-20, maxIter = 2)


##################################################
## Distance matrix for places
## install.packages("maps")
siteNstate <- paste(workingDf$Site, workingDf$State)
unique(siteNstate)
write.csv(sites, file = "sites.csv")





##################################################
## try phangorn
library(ape)
library(phangorn)
library(multicore)
primates <- read.phyDat("primates.dna", format="phylip", type="DNA")
dm <- dist.dna(as.DNAbin(primates))
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")

fit <- pml(treeNJ, data=primates)
fit$weight
trCl <- list()
class(trCl) <- "phylo"
trCl$Nnode <- 3
tc$edge <- matrix(c(5,1,5,2,6,3,6,4,7,5,7,6),ncol=2 , byrow=TRUE)
trCl$tip.label <- c("A","B","C","D", "E", "F")
plot(trCl, "cladogram")


Ntip <- length(treeNJ$tip.label)
Nedge <- dim(treeNJ$edge)[1]
Nnode <- treeNJ$Nnode
ROOT <- Ntip + 1
yy <- numeric(Ntip + Nnode)
pc <- pml(tree = tc, data = sc)



library(phangorn)
library(abind)
library(inline)
library(RcppArmadillo)
library(nnls)


postRcpp <- cxxfunction(signature(data ="list", P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines( "likelihood.cpp"), collapse = "\n" ))




set.seed(42)
tc <- read.tree(text = "((a,b,c),(d,e,f));")
tc$edge.length <- c(0.3,0.1,0.1,0.1,0.3,0.1,0.1,0.1)
sc <- simSeq(tc, l=30)
a <- pcl(tc, sc)
P <- a$P[1,1][[1]]
contrast <- attr(sc, "contrast")
nrs <- as.integer(length(sc$a))
ncs <- as.integer(attr(sc,"nc"))
ecps <- matrix(c(0.6,0.6,0.4,0.4,0.4,0.4,0.6,0.6), 4, 2)
initP <- matrix(c(rep(0.55,3),rep(0.45,3), rep(0.45,3), rep(0.55,3)), 6,2)


b <- postRcpp(data=sc[tc$tip.label], P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=18, bfs=c(0.25,0.25,0.25,0.25), ecps=initP)


sDNAc <- as.DNAbin(sc)
distc <- as.matrix(dist.dna(sDNAb))
y <- distb[lower.tri(distb)]
dataTree <- designObsClade(distb)

## try this
debug(mtCluster)
mtCluster(sc, P = P, distMat = distb, M = 2, initP=initP, maxIter = 50, tol = 1e-10)


plot(tc, "cladogram")
ecps <- matrix(c(0.6,0.6,0.4,0.4,0.4,0.4,0.6,0.6), 4, 2)
system.time(temp <- seqEM(tree=tc, seqData=sc, P=P, init=ecps2, iter=50, method="cpp"))




r <- list()
system.time(r <- simseqEM(tc=tc[[2]],P=P, iter=100))


tc <- list()
tc[[1]] <- read.tree(text = "((a,b,c,d,e),(f,g,h,i));")
tc[[2]] <- read.tree(text = "((a,b,c),(d,e,f),(g,h,i));")
set.seed(42)


nCorrect(tc[[2]], r, 0.8)

set.seed(42)
n75 <- n90 <- n95 <- nC <- matrix(,6,2)
edges <- seq(0.05,0.3,by=0.05)
for(i in 1:length(edges)){
    for(j in 1:2){
        r <- simseqEM(tc=tc[[j]],P=P, edgeLength = edges[i], iter=500)
        nC[i,j] <- nCorrect(tc[[j]], r, 0)
        n75[i,j] <- nCorrect(tc[[j]], r, 0.75)
        n90[i,j] <- nCorrect(tc[[j]], r, 0.9)
        n95[i,j] <- nCorrect(tc[[j]], r, 0.95)
    }
}


sum(sc[[1]]==sc[[2]])
sum(sc[[1]]==sc[[3]])
sum(sc[[1]]==sc[[4]])
sum(sc[[2]]==sc[[4]])
sum(sc[[2]]==sc[[3]])
sum(sc[[4]]==sc[[3]])


nrs <- as.integer(length(sc$a))
ncs <- as.integer(attr(sc,"nc"))
contrast <- attr(sc, "contrast")


a <- pcl(tc, sc)
P <- a$P[1,1][[1]]
contrast <- a$contrast
b <- contrast %*% P[1,1][[1]]
d <- list(rep(1,5),rep(2,5),rep(3,5))

temp <- read.tree(text = "((a,b),(c, d));")
