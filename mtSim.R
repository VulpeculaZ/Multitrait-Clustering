library(phangorn)
## library(abind)
library(inline)
library(RcppArmadillo)
## Cannot work in R2.14
## library(openNLP)

tb <- read.tree(text = "((a,b,c),(d,e,f),(g,h));")
tb <- reorderPruning(tb)
tb$edge.length <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2)
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





## Examples:

a <- simClSeq(2, 6, group = c(3,3), outLen = 0.1, inLen = 0.2)
a <- simDistMat(2, 6, group = c(3,3), outLen = 0.1, inLen = 0.2, sig = 0.01)


##################################################
## Simulation: Comparing using one trait alone vs. using two
##################################################

aSeq <- simClSeq(2, 10, group = c(5,5), outLen = 0.1, inLen = 0.2, returnTree =TRUE)
tc <- aSeq$tree
a <- pcl(aSeq$tree, aSeq$seq)
P <- a$P[1,1][[1]]


set.seed(42)
simResult <- list()
convResult <- list()
for(i in 1:20){
    initP <- randomPMat(pDiff = 0.05, 2, 10)
    aSeq <- simClSeq(2, 10, group = c(5,5), outLen = 0.3, inLen = 0.1, l=200)
    aDist <- simDistMat(2, 10, group = c(5,5), outLen = 0.3, inLen = 0.1, sig = 0.5)
    tryCatch({tempResult <- mtCluster(seqData = aSeq, P = P, distMat = aDist, M=2, initP = initP, maxIter = 100, tol = 1e-4); simResult[[i]] <- tempResult$Ez; convResult[[i]] <- tempResult$convergence},
             error =function(e)
         {a <<- list(initP=initP, seq = aSeq, dist =aDist)
          print(i)}
             )
}


nCorrect(tc = tc, r=simResult)
monotoneEM(convResult)


## is DNA monotone in EM???

set.seed(42)
simResult <- list()
convResult <- list()
for(i in 1:100){
    initP <- randomPMat(pDiff = 0.05, 2, 10)
    aSeq <- simClSeq(2, 10, group = c(5,5), outLen = 0.3, inLen = 0.1, l=200)
    tempResult <- DNACluster(seqData = aSeq, P = P, M=2, initP = initP, maxIter = 100, tol = 1e-4)
    simResult[[i]] <- tempResult$Ez;
    convResult[[i]] <- tempResult$convergence
}

## is DistMat monotone in EM???
set.seed(42)
simResult <- list()
convResult <- list()

for(i in 1:20){
    initP <- randomPMat(pDiff = 0.05, 2, 10)
    aDist <- simDistMat(2, 10, group = c(5,5), outLen = 0.1, inLen = 0.2, sig = 0.1)
    tryCatch({tempResult <- DistCluster(distMat = aDist, M=2, initP = initP, maxIter = 100, tol = 1e-4); simResult[[i]] <- tempResult$Ez.old; convResult[[i]] <- tempResult$convergence},
             error =function(e)
         {a <<- list(initP=initP, dist =aDist)
          print(i)}
             )
}

nCorrect(tc = tc, r=simResult)
monotoneEM(convResult)

initP <- randomPMat(pDiff = 0.05, 2, 10)
aDist <- simDistMat(2, 10, group = c(5,5), outLen = 0.2, inLen = 0.3, sig = 0.1)

initP <- matrix(c(rep(0.1,5), rep(0.9,10), rep(0.1,5)), 10,2)

set.seed(42)
tempResult <- DistCluster(distMat = aDist, M=2, initP = initP, maxIter = 4, tol = 1e-4)

## New likelihood for hosts!!!


distMat <- matrix(c(0,1,2,2,
                    1,0,2,2,
                    2,2,0,1,
                    2,2,1,0),
                  4,4)

set.seed(42)
eObs <- c(sample(c(1,2), 6, replace = TRUE), sample(c(3,4),6,replace =TRUE))
error <- rgeom(12, prob=0.6)
error[error > 2] <- 2
obs <- apply(cbind(distMat[eObs,], error), 1, FUN = rObs)
initP <- randomPMat(pDiff = 0.05, 2, 12)
initQ <- randomPMat(pDiff = 0.05, 2, 4)
hostCluster(obs, distMat, initP, maxIter = 200)

## Test the new host likelihood with sequences.
distMat <- matrix(c(0,1,2,2,
                    1,0,2,2,
                    2,2,0,1,
                    2,2,1,0),
                  4,4)

set.seed(42)
eObsHost <- c(sample(c(1,2), 6, replace = TRUE), sample(c(3,4),6,replace =TRUE))
errorHost <- rgeom(12, prob=0.6)
errorHost[errorHost > 2] <- 2
obsHost <- apply(cbind(distMat[eObsHost,], errorHost), 1, FUN = rObs)
initP <- randomPMat(pDiff = 0.05, 2, 12)


## Simulation for multiple trait:
## One simple example:
obsSeq <- simClSeq(2, 12, group = c(6,6), outLen = 0.3, inLen = 0.1, l=200)
mtCluster(obsSeq, obsHost=obsHost, distMat=distMat, initP = initP)
##


##################################################
## Simulation for cross-validation
## A simple CV:
distMat <- matrix(c(0,1,2,3,3,
                    1,0,2,3,3,
                    2,2,0,3,3,
                    3,3,3,0,1,
                    3,3,3,1,0),
                  5,5)


set.seed(42)
eObsHost <- c(sample(c(1,2,3), 15, replace = TRUE), sample(c(4,5),15,replace =TRUE))
errorHost <- rgeom(30, prob = 0.6)
errorHost[errorHost > 3] <- 3
obsHost <- apply(cbind(distMat[eObsHost,], errorHost), 1, FUN = rObs)
obsSeq <- simClSeq(2, 30, group = c(15,15), outLen = 0.3, inLen = 0.1, l=200)
initP <- randomPMat(pDiff = 0.05, 1, 30)


cvDataList <- list(obsSeq=obsSeq, P=P, obsHost=obsHost, distMat = distMat)
cvll <- cvCluster(cvDataList, 30, maxCluster = 3, mCV=100, beta=0.5)


##################################################
## CV: 2 clusters are prefered:
distMat <- matrix(c(0,1,2,2,
                    1,0,2,2,
                    2,2,0,1,
                    2,2,1,0),
                  4,4)

set.seed(42)
cvll <- list()
for(i in 1:100){
    eObsHost <- c(sample(c(1,2), 15, replace = TRUE), sample(c(3,4),15,replace =TRUE))
    errorHost <- rgeom(30, prob = 0.6)
    errorHost[errorHost > 2] <- 2
    obsHost <- apply(cbind(distMat[eObsHost,], errorHost), 1, FUN = rObs)
    obsSeq <- simClSeq(2, 30, group = c(15,15), outLen = 0.3, inLen = 0.1, l=200)
    initP <- randomPMat(pDiff = 0.05, 1, 30)
    cvDataList <- list(obsSeq=obsSeq, P=P, obsHost=obsHost, distMat = distMat)
    cvll[[i]] <- cvCluster(cvDataList, 30, maxCluster = 3, mCV=100, beta=0.5)
}

cvSim <- rep(0,3)
for(i in 1:100){
    cvSim <- cvSim + (cvll[[i]] == max(cvll[[i]]))
}

