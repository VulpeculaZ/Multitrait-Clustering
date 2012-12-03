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
# Simulation: Comparing using one trait alone vs. using two
#################################################

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


## Testing clustering for DNA sequences.

set.seed(42)
simResult <- list()
convResult <- list()

nSeq <- c(25, 50, 100)
edgeLength <- c(0.2, 0.5, 1)
nCluster <- c(2,3,4)

for(iSeq in nSeq){
    for(iedge in edgeLength)
        for(iCluster in nCluster){
            for(i in 1:100){
                initP <- randomPMat(pDiff = 0.05, iCluster, iSeq)
                iGroup <- table(sample(1:iCluster, iSeq - 5 * iCluster, replace =TRUE))
                aSeq <- simClSeq(iCluster, iSeq, group = c(5,5), outLen = 0.5, inLen = 0.5, l=20)
                tempResult <- DNACluster(seqData = aSeq, initP = initP, maxIter = 100, tol = 1e-4, P=P, method = "EMmn")
}

        }
    }
}

## is DNA monotone in EM???

for(i in 1:100){
    initP <- randomPMat(pDiff = 0.05, 2, 10)
    aSeq <- simClSeq(2, 10, group = c(5,5), outLen = 0.5, inLen = 0.5, l=20)
    tempResult <- DNACluster(seqData = aSeq, initP = initP, maxIter = 100, tol = 1e-4, P=P, method = "EMmn")
    simResult[[i]] <- tempResult$Ez;
    convResult[[i]] <- tempResult$convergence
}

monotoneEM(convResult)

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
hostCluster(obs, distMat, initP, maxIter = 200, model="normal")

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

##################################################
## two moon
##################################################
library(spa)

set.seed(100)
twoMoons <- spa.sim(type="moon")
plot(twoMoons$x1, twoMoons$x2)
distMat <- sqrt(outer(twoMoons$x1, twoMoons$x1, FUN = "-") ^ 2 + outer(twoMoons$x2, twoMoons$x2, FUN = "-") ^ 2)


twoC <- list()
twoC$x1 <- c(rnorm(103), rnorm(103) + 3)
twoC$x2 <- c(rnorm(103), rnorm(103) + 3)
distMat2c <- sqrt(outer(twoC$x1, twoC$x1, FUN = "-") ^ 2 + outer(twoC$x2, twoC$x2, FUN = "-") ^ 2)
plot(twoMoons$x1, twoMoons$x2)


## initP <- matrix(c(rep(0.7,103), rep(0.3,206), rep(0.7,103)),206,2)
initP <- randomPMat(pDiff = 0.1, 2, 206)
mmTm2 <- mmEm(distMat, initP =initP, tol=1e-08, maxIter = 1)
notCorrect <- as.logical(abs((mmTm$Ez < 0.5) [,1]  - c(rep(1,103), rep(0,103))))
points(twoMoons$x1[notCorrectNormal], twoMoons$x2[notCorrectNormal], col = 2)
notCorrectNormal <- as.logical(4* twoMoons$x1 -6 - twoMoons$x2 > 0)
notCorrectNormal[104:206] <- FALSE

##' ##' <description>
##'
##' <details>
##' @title
##' @param y1 Data
##' @param y2 Data
##' @param x1 Prediction
##' @param x2 Prediction
##' @param t
##' @param theta
##' @return
##' @author Ziqian Zhou
mmEmPred <- function(y1, y2, x1, x2, t, theta){
    M <- dim(t)[2]
    n <- dim(t)[1]
    distXY <- sqrt(outer(x1, y1, FUN = "-") ^ 2 + outer(x2, y2, FUN = "-") ^ 2)
    distYY <- sqrt(outer(y1, y1, FUN = "-") ^ 2 + outer(y2, y2, FUN = "-") ^ 2)
    predP <- matrix(, nrow = length(x1), ncol=M)
    dMat <- dnorm(x = distXY, mean =0, sd=theta)
    pMat <- dnorm(x = distYY, mean = 0, sd = theta)

    for(i in 1:M){
        pMatt <- sweep(pMat, 2, t[,i], FUN = "^")
        qi <- apply(pMatt,1, prod)
        qi <- qi / sum(qi)
        predP[,i] <- rowSums(sweep(dMat, 2, qi, FUN = "*"))
    }
    predP <- sweep(predP,1, rowSums(predP), "/")
    predP
}

xPred <- list(x1 =c(2,4,3,3.08), x2 = c(2,-4,-2,1.345))

x1Pred <-  seq(0,6,len = 200)
x2Pred <- seq(-4,2,len = 200)
xpredList <- list(x1 = rep(x1Pred, each=200) , x2 = rep(x2Pred, 200))
dContour <- mmEmPred(twoMoons$x1, twoMoons$x2,xpredList$x1, xpredList$x2, mmTm$Ez, 0.5)
mContour1 <- matrix(dContour[,1], 200,200)
mContour2 <- matrix(dContour[,2], 200,200)

mContour1 <- mContour1-0.99
mContour1[mContour1<0] <-  0
contour(mContour1)

lines(x =c(1,2), y = c(-2,2))
twoMoons$x1
