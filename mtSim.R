library(phangorn)
## library(abind)
library(inline)
library(RcppArmadillo)
## Cannot work in R2.14
## library(openNLP)
library(nnls)

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

##' <description>
##'Create random initial probablity matrix
##' <details>
##' @title Random initial probablity matrix
##' @param pDiff Difference in probality.
##' @param m Number of clusters
##' @param n Number of observations
##' @return A matrix.
##' @author Ziqian Zhou
randomPMat <- function(pDiff, m, n){
    if(abs(pDiff) >1 ) stop("abs(pDiff) must be smaller than 1" )
    p <- runif(m*n, (1-pDiff), (1+pDiff))
    pMat <- matrix(p, n, m)
    pMat <- pMat / rowSums(pMat)
}

##' <description>
##' Simulate sequences base on cluster grouping.
##' <details>
##' @title Make Sequecnes from Clusters
##' @param M Number of groups
##' @param n Number of observations
##' @param group The size of groups. Should be a vector of integer of
##' size M.
##' @param outLen The length of outter branches. Should  be a vector
##' of size n or 1.
##' @param inLen The length of internal branches. Should be a vector
##' of size M or 1.
##' @param l The length of each sequence.
##' @param returnTree Whether or not return the tree.
##' @return A group of n sequences of the class phyDat.
##' @author Ziqian Zhou
simClSeq <- function(M, n, group, outLen, inLen, l = 500, returnTree = FALSE){
    rTree <- list()
    rTree$edge <- matrix(, nrow =n+M, ncol =2)
    rTree$edge[1:n, 1] <- rep((n+2):(n+M+1), group)
    pos <- 1
    for(i in 1:M){
        rTree$edge[pos:(pos+group[i]-1), 2] <- (pos+group[i]-1):pos
        pos <- pos + group[i]
    }
    rTree$edge[(n+1) : (n+M), 1] <- n+1
    rTree$edge[(n+1) : (n+M), 2] <- (n+M+1) : (n+2)
    rTree$edge.length <- NA
    rTree$edge.length[1:n] <- outLen
    rTree$edge.length[(n+1):(n+M)] <- inLen
    rTree$Nnode <- as.integer(M + 1)
    rTree$tip.label <- as.character(1:n)
    attr(rTree, "class") = "phylo"
    attr(rTree, "order") = "pruningwise"
    if(!isTRUE(returnTree)){
        return(simSeq(rTree, l = l))
    }
    else{
        return(list(seq=simSeq(rTree, l=l), tree = rTree))
    }
}




##' <description>
##' Simulate distance matrix base on cluster grouping.
##' <details>
##' @title Make distance matrix from Clusters
##' @param M Number of groups
##' @param n Number of observations
##' @param group The size of groups. Should be a vector of integer of
##' size M.
##' @param outLen The length of outter branches. Should  be a vector
##' of size n or 1.
##' @param inLen The length of internal branches. Should be a vector
##' of size M or 1.
##' @param sig The standard deviation multiplier.
##' @return A distance matrix.
##' @author Ziqian Zhou
simDistMat <- function(M, n, group, outLen, inLen, sig){
    if(length(inLen)==1) inLen <- rep(inLen, M)
    if(length(outLen)==1) outLen <- rep(outLen, n)
    matTrue <- matrix(0, n, n)
    ipos <- 1
    jpos <- 1 + group[1]
    for(i in 1:(M-1)){
        for(j in (i+1):M){
            matTrue[ipos:(ipos+group[i]-1), jpos:(jpos+group[j]-1)] <- inLen[i] + inLen[j]
            ipos <- ipos + group[i]
            jpos <- jpos + group[j]
        }
    }
    matTrue[lower.tri(matTrue)] <- t(matTrue)[lower.tri(matTrue)]
    matTrue <- matTrue + outer(outLen, outLen, FUN = "+")
    error <- matrix(rnorm(n*n, sd = sig * sqrt(matTrue)), n, n)
    error[lower.tri(error)] <- t(error)[lower.tri(error)]
    distMat <- matTrue + error
    diag(distMat) <- 0
    pmax(distMat, 0)
}


##' <description>
##' Test the number of correctly identified clusters.
##' <details>
##' @title nCorrect
##' @param tc The true tree.
##' @param r A list of matrices of classification probabilities.
##' @param p Cut off probabiltiy for classification. If p==0, use the cluster with highest probability as the classification result.
##' @return The ratio of correctly identified clusters.
##' @author Ziqian Zhou
nCorrect <- function(tc, r, p=0){
    nTips <- length(tc$tip.label)
    inner <- (nTips+2):(nTips+tc$Nnode)
    groups <- list()
    ratioC <- 0
    for(i in 1:length(inner)){
        groups[[i]] <- sort(as.integer(tc$edge[which(tc$edge[,1] == inner[i]), 2]))
    }
    for(i in 1:length(r)){
        idr <- list()
        maxr <- apply(r[[i]], 1, function(x) which(x==max(x)))
        for(j in 1:length(inner)){
            idr[[j]] <- sort(as.integer(which(maxr == j)))
            for(k in 1:length(inner)){
                if(p==0) ratioC <- ratioC + identical(idr[[j]], groups[[k]])
                else ratioC <- ratioC + as.integer(identical(idr[[j]], groups[[k]]) && min(r[[i]][idr[[j]],j]) > p)
            }
        }
    }
    return(ratioC / length(inner) / length(r))
}



## Examples:

a <- simClSeq(2, 6, group = c(3,3), outLen = 0.1, inLen = 0.2)
a <- simDistMat(2, 6, group = c(3,3), outLen = 0.1, inLen = 0.2, sig = 0.01)


##################################################
## Simulation: Comparing using one trait alone vs. using two
##################################################

postRcpp <- cxxfunction(signature(data ="list", P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/likelihood.cpp"), collapse = "\n" ))
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


monotoneEM <- function(x){
    mono <- 0
    for(i in 1:length(x)){
        if(any(x[[i]][-1] - x[[i]][-length(x[[i]])] < 0)){
            mono <- mono + 1
            print(i)
        }
    }
    mono
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
aDist <- simDistMat(2, 10, group = c(5,5), outLen = 0.1, inLen = 0.3, sig = 0.1)

initP <- matrix(c(rep(0.1,5), rep(0.9,10), rep(0.1,5)), 10,2)

tempResult <- DistCluster(distMat = aDist, M=2, initP = initP, maxIter = 100, tol = 1e-4)

