library(phangorn)
## library(abind)
library(inline)
library(RcppArmadillo)
## library(openNLP)
library(combinat)
library(parallel)
library(cluster)
library(nnet)

tb <- read.tree(text = "((a,b,c),(d,e,f),(g,h));")
tb <-  phangorn:::reorderPruning(tb)
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

aSeq <- simClSeq(2, 6, group = c(3,3), outLen = 1, inLen = 1, returnTree =TRUE)
tc <- aSeq$tree
tc$tip.label <- c("a","b", "c", "d", "e", "f")
plot(tc)
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
dimnames = c("metho", "nSeq", "nCluster", "edgeLength")
ratioCorrect <- array(dim = c(2,3,3,4))
meanTheta <- sdTheta<- array(dim = c(2,3,3,4))
thetaAll <- list()

convResult <- list()
nSeq <- c(24, 48, 96)
edgeLength <- c(0.2, 0.5, 1, 2)
nCluster <- c(2,3,4)
edgeRatio <- c(1,2)
j <- 1

for(iRatio in edgeRatio){
for(iSeq in nSeq){
    for(iCluster in nCluster){
        for(iedge in edgeLength){
            iGroup <- rep(iSeq/iCluster, iCluster)
            resultSeqEz <- list()
            resultMNEz <- list()
            resultSeqTheta <- list()
            resultMNTheta <- list()
            for(i in 1:100){
                tempSeq <- simClSeq(iCluster, iSeq, group = iGroup, outLen = iedge, inLen = iedge * iRatio, l=100)
                tempDNA <- as.DNAbin(tempSeq)
                distMat <- as.matrix(dist.dna(tempDNA, model = "raw"))
                initP <- EMinit(distMat, iCluster, diffP = 0.1)
                tempSeqEM <- DNACluster(seqData = tempSeq, initP = initP, maxIter = 100, tol = 1e-5)
                resultSeqEz[[i]] <- tempSeqEM$Ez
                resultSeqTheta[[i]] <- tempSeqEM$theta
                if(any(is.na(tempSeqEM$theta))){
                    badSeq <- tempSeq
                    badInitP <- initP
                    break()
                }
            }
            thetaAll[[j]] <- unlist(resultSeqTheta)
            j <- j + 1
            ratioCorrect[iRatio, which(iSeq==nSeq), which(iCluster==nCluster), which(iedge==edgeLength)] <- nCorrect(iGroup, resultSeqEz)
            meanTheta[iRatio, which(iSeq==nSeq), which(iCluster==nCluster), which(iedge==edgeLength)] <- mean(unlist(resultSeqTheta),na.rm = TRUE)
            sdTheta[iRatio, which(iSeq==nSeq), which(iCluster==nCluster), which(iedge==edgeLength)] <- sd(unlist(resultSeqTheta), na.rm = TRUE)
        }
    }
}
}


## Testing the distribution of theta.

set.seed(42)
dimnames = c("nCluster", "edgeLength")
sdTheta <- 0
thetaAll <- list()
edgeLength <- c(0.2, 0.5, 1)
iRatio <- 2
iSeq <- 96
j <- 1

for(iCluster in nCluster){
    for(iedge in edgeLength){
        iGroup <- rep(iSeq/iCluster, iCluster)
        for(i in 1:100){
            tempSeq <- simClSeq(iCluster, iSeq, group = iGroup, outLen = iedge, inLen = iedge * iRatio, l=100)
            tempDNA <- as.DNAbin(tempSeq)
            distMat <- as.matrix(dist.dna(tempDNA, model = "raw"))
            initP <- EMinit(distMat, iCluster, diffP = 0.1)
            tempSeqEM <- DNACluster(seqData = tempSeq, initP = initP, maxIter = 100, tol = 1e-5)
            resultSeqEz[[i]] <- tempSeqEM$Ez
            resultSeqTheta[[i]] <- tempSeqEM$theta
        }
        thetaAll[[j]] <- unlist(resultSeqTheta)
        sdTheta[j] <- sd(unlist(resultSeqTheta))
        j <- j + 1
    }
}

thetaStd <- list()

par(mfrow=c(3,3))

for(j in 1:9){
    thetaStd[[j]] <- (thetaAll[[j]] - mean(thetaAll[[j]])) / sdTheta[j]
    qqnorm(thetaStd[[j]], main=NULL)
    qqline(thetaStd[[j]])
}

qqnorm
## boot simulation


tempSeqEM <- DNACluster(seqData = badSeq, initP = initP, maxIter = 100, tol = 1e-4,  method = "EMseq")


## is DNA monotone in EM???

for(i in 1:100){
    initP <- randomPMat(pDiff = 0.05, 2, 10)
    aSeq <- simClSeq(2, 10, group = c(5,5), outLen = 0.5, inLen = 0.5, l=50)
    tempResult <- DNACluster(seqData = aSeq, initP = initP, maxIter = 100, tol = 1e-4,  method = "EMseq")
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
initP <- matrix(c(rep(0.55, 6), rep(0.45, 12), rep(0.55, 6)),12,2)
hostCluster(obs, distMat, initP, maxIter = 200, model="geom")

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
## Test for speed
obsSeq <- simClSeq(2, 500, group = c(250,250), outLen = 0.1, inLen = 0.3, l=500)
initP <- randomPMat(pDiff = 0.05, 2, 500)
system.time(tempSeqEM <- DNACluster(seqData = obsSeq, initP = initP, maxIter = 100, tol = 1e-5,  method = "EMseq"))
system.time(tempSeqEMRcpp <- DNACluster(seqData = obsSeq, initP = initP, maxIter = 100, tol = 1e-5,  method = "cpp"))


##################################################
## Initalization using k-Medoids
obsSeq <- simClSeq(2, 1000, group = c(500, 500), outLen = 0.1, inLen = 0.3, l=500)
seqDNA <- as.DNAbin(obsSeq)
distMat <- as.matrix(dist.dna(seqDNA))

debug(EMinit)
initP <- EMinit(distMat, 2)

##################################################
## Test bootstrap
obsSeq <- simClSeq(2, 500, group = c(250,250), outLen = 0.1, inLen = 0.3, l=500)
initP <- randomPMat(pDiff = 0.05, 2, 500)
tempSeqEM <- DNACluster(seqData = obsSeq, initP = initP, maxIter = 100, tol = 1e-5,  method = "cpp")
theta <- tempSeqEM$theta
y <- tempSeqEM$y
Ez <- tempSeqEM$Ez
obs <- obsSeq

##################################################
## bootstrap simulation
set.seed(42)
initP <- matrix(c(rep(0.6,100), rep(0.4, 200), rep(0.6,100)),200,2 )
alphaMin <- alphaQuad <- 0

for(i in 1:10){
    obsSeq <- simClSeq(2, 400, group = c(200,200), outLen = 0.1, inLen = 0.3, l=500)
    train <- obsSeq[c(1:100, 301:400)]
    test <- obsSeq[101:300]
    attr(train, "contrast") <- attr(obsSeq, "contrast")
    attr(train, "nc") <- attr(obsSeq, "nc")
    tempSeqEM <- DNACluster(seqData = train, initP = initP, maxIter = 100, tol = 1e-5,  method = "cpp")
    theta <- tempSeqEM$theta
    p <- tempSeqEM$p
    y <- tempSeqEM$y
    Ez <- tempSeqEM$Ez
    bTheta <- bootTheta(theta=theta, y=y, p=p,obs = train, Ez = Ez ,n=10000, percentile = 0.95)
    cutQuad <- bTheta$qQuad
    varScore <- bTheta$varScore
    boot <- bootX(theta=theta, y=y, p=p, x=test, varScore = varScore)
    ## sum(boot$bootMin > cutMin)
    sum(boot$bootQuad > cutQuad)
    ## alphaMin <- alphaMin + sum(boot$bootMin > cutMin)
    alphaQuad <- alphaQuad + sum(boot$bootQuad > cutQuad)
}

set.seed(42)
convResult <- list()
nSeq <- c(24, 48, 96)
edgeLength <- c(0.2, 0.5, 1)
nCluster <- c(2,3,4)
edgeRatio <- c(1,2)
initP <- matrix(c(rep(0.6,100), rep(0.4, 200), rep(0.6,100)),200,2 )
alpha <- alphaN <- matrix(0,2,3)
testPower <- testPowerN <- matrix(0,2,3)
cutChi <- qchisq(0.95, 1) / 200

for(iRatio in edgeRatio){
    for(iedge in edgeLength){
        for(i in 1:100){
            obsSeq <- simClSeq(3, n = 600, group = c(200,200,200), outLen = iedge, inLen = iedge * iRatio, l=100, ancestral = TRUE)
            train <- obsSeq[c(1:100, 301:400)]
            test <- obsSeq[101:300]
            out <- obsSeq[401:600]
            ans <- obsSeq[601:604]
            attr(train, "contrast") <- attr(obsSeq, "contrast")
            attr(train, "nc") <- attr(obsSeq, "nc")
            tempSeqEM <- DNACluster(seqData = train, initP = initP, maxIter = 100, tol = 1e-8,  method = "cpp")
            theta <- tempSeqEM$theta
            Ez <- tempSeqEM$Ez
            p <- tempSeqEM$p
            y <- tempSeqEM$y
            bTheta <- bootTheta(theta=theta, y=y, p=p,obs=train, Ez = Ez, n=10000, percentile = 0.95)
            ## bootAns <- bootX(theta=theta, y=y, p=p, x=ans, varScore = varScore)
            cutQuad <- bTheta$qQuad
            varScore <- bTheta$varScore
            boot <- bootX(theta=theta, y=y, p=p, x=test, varScore = varScore)
            bootOut <- bootX(theta=theta, y=y, p=p, x=out, varScore = varScore)
            ## bootArr[iRatio,  which(iedge==edgeLength), ] <- boot
            alpha[iRatio,  which(iedge==edgeLength)] <- alpha[iRatio, which(iedge==edgeLength)] + sum(boot$bootQuad > cutQuad)
            testPower[iRatio, which(iedge==edgeLength)] <- testPower[iRatio,  which(iedge==edgeLength)] +  sum(bootOut$bootQuad > cutQuad)
            alphaN[iRatio,  which(iedge==edgeLength)] <- alphaN[iRatio, which(iedge==edgeLength)] + sum(boot$bootQuad > cutChi)
            testPowerN[iRatio, which(iedge==edgeLength)] <- testPowerN[iRatio,  which(iedge==edgeLength)] +  sum(bootOut$bootQuad > cutChi)
        }
    }
}


##################################################
## CV for DNA using EM
initPList <- list()


cvDataList <- list()
cvll <- list()
initPList <- list()
iter <- c(2:4)

set.seed(42)

for(i in 1:1000){
    cvDataList$obsSeq <- simClSeq(3, 50, group = c(20,15,15), outLen = 0.3, inLen = 0.1, l=200)
    for(j in iter){
        initPList[[j]] <- randomPMat(pDiff = 0.05, j, 50)
    }
    cvll[[i]] <- cvCluster(cvDataList, 50, maxCluster = c(2:4), mCV=100, beta=0.5, initPList = initPList)
}

cvResult <- rep(0,3)
for(i in 1:1000){
    colSumsCv <- colSums(cvll[[i]])
    cvResult <- cvResult + (colSumsCv == max(colSumsCv))
}

cvResult2 <- cvResult1 <- rep(0,3)
for(i in 1:1000){
    maxCv <- rowSums(apply(cvll[[i]], 1, function(x) x == max(x)))
    cvResult1 <- cvResult1 + maxCv
    cvResult2 <- cvResult2 + (maxCv == max(maxCv))
}

nCluster <- 3
nSeq <- seq(40, 80, by = 10)
nLen <- c(250, 500, 1000)
diffResult <- matrix(0,length(nSeq), length(ratio))
nPara <- outer(nSeq, ratio)


set.seed(42)
for(i in nSeq){
    for(j in ratio){
        for(k in 1:1000){
            tempGroup <- i * c(0.3, 0.3, 0.4)
            tempInitP <- randomPMat(0.05, nCluster, i)
            allSeq <- simClSeq(nCluster, i, group = tempGroup, outLen = 0.5, inLen = 0.5, l= floor(i * j), ancestral = TRUE)
            anceSeq <- allSeq[(i+2):(i+nCluster+1)]
            anceSeq <- matrix(unlist(anceSeq), ncol=nCluster)
            tempSeq <- allSeq[1:i]
            attr(tempSeq, "contrast") <- attr(allSeq, "contrast")
            attr(tempSeq, "nc") <- attr(allSeq, "nc")
            tempResult <- DNACluster(tempSeq, initP = tempInitP, maxIter = 100, tol = 1e-5,  method = "EMseq")
            diffResult[which(i == nSeq), which(j == ratio)] <- diffResult[which(i == nSeq), which(j == ratio)] + mDiff(anceSeq, tempResult$y)
        }
    }
}

nSeq <- seq(40, 80, by = 10)
nLen <- c(250, 500, 1000, 2000)
testForRate <- matrix(0, length(nSeq), length(nLen))

for(j in 1:length(nLen)){
    iLen <- nLen[j]
    for(i in 1:length(nSeq)){
    iSeq <- nSeq[i]
    for(k in 1:100){
        iGroup <- c(0.5, 0.5) * iSeq
        initP <- matrix(c(rep(0.55, iSeq /2), rep(0.45, iSeq), rep(0.55, iSeq/2)), iSeq, 2)
        allSeq <- simClSeq(2, iSeq, iGroup, outLen= 0.5, inLen = 0.5, l=iLen, ancestral = TRUE)
        anceSeq <- allSeq[(iSeq+2):(iSeq+3)]
        anceSeq <- matrix(unlist(anceSeq), ncol=2)
        tempSeq <- allSeq[1:iSeq]
        attr(tempSeq, "contrast") <- attr(allSeq, "contrast")
        attr(tempSeq, "nc") <- attr(allSeq, "nc")
        tempResult <- DNACluster(tempSeq, initP = initP, maxIter = 100, tol = 1e-5,  method = "EMseqRcpp")
        testForRate[i,j]<-  testForRate[i,j] + mDiff(anceSeq, tempResult$y)
    }
}
}

diffResult / 1000 /3

mDiff <- function(ance, result){
    m <- dim(ance)[2]
    nm <- m * dim(ance)[1]
    allPerm <- permn(x =1:m)
    nm - max(unlist(lapply(allPerm, function(x) sum(ance[,x] == result))))
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

