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
simClSeq <- function(M, n, group, outLen, inLen, l = 500, returnSeq = TRUE, returnTree = FALSE, ancestral = FALSE){
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
    if(isTRUE(returnSeq)){
        if(!isTRUE(returnTree)){
            return(simSeq(rTree, l = l, ancestral = ancestral))
        }
        else{
            return(list(seq=simSeq(rTree, l=l,  ancestral = ancestral), tree = rTree))
        }
    }
    else{
        return(rTree)
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
## nCorrect <- function(tc , r, p=0){
##     nTips <- length(tc$tip.label)
##     inner <- (nTips+2):(nTips+tc$Nnode)
##     groups <- list()
##     ratioC <- 0
##     for(i in 1:length(inner)){
##         groups[[i]] <- sort(as.integer(tc$edge[which(tc$edge[,1] == inner[i]), 2]))
##     }
##     for(i in 1:length(r)){
##         idr <- list()
##         maxr <- apply(r[[i]], 1, function(x) which(x==max(x)))
##         for(j in 1:length(inner)){
##             idr[[j]] <- sort(as.integer(which(maxr == j)))
##             for(k in 1:length(inner)){
##                 if(p==0) ratioC <- ratioC + identical(idr[[j]], groups[[k]])
##                 else ratioC <- ratioC + as.integer(identical(idr[[j]], groups[[k]]) && min(r[[i]][idr[[j]],j]) > p)
##             }
##         }
##     }
##     return(ratioC / length(inner) / length(r))
## }


##' <description>
##' Test the number of correctly identified clusters.
##' <details>
##' @title nCorrect
##' @param tc The true tree.
##' @param r A list of matrices of classification probabilities.
##' @param p Cut off probabiltiy for classification. If p==0, use the cluster with highest probability as the classification result.
##' @return The ratio of correctly identified clusters.
##' @author Ziqian Zhou
nCorrect <- function(grouping, r, p=0){
    n <- length(r)
    nr <- dim(r[[1]])[1]
    ndiv <- 0
    for(j in 1:length(r)){
        rGrouping <- max.col(r[[j]])
        tGrouping <- table(rGrouping)
        cumsumGrouping <- cumsum(grouping)
        for(i in 1:length(grouping)){
            if(i == 1){
                iNGroup <- max(table(rGrouping[1:cumsumGrouping[i]]))
                iGroup <- which.max(table(rGrouping[1:cumsumGrouping[i]]))
                ndiv <- ndiv + abs(iNGroup - grouping[i]) + abs(tGrouping[iGroup] - iNGroup)
            }
            else{
                iNGroup <- max(table(rGrouping[cumsumGrouping[i-1]:cumsumGrouping[i]]))
                iGroup <- which.max(table(rGrouping[cumsumGrouping[i-1]:cumsumGrouping[i]]))
                ndiv <- ndiv + abs(iNGroup - grouping[i]) + abs(tGrouping[iGroup] - iNGroup)
            }
        }
    }
    return(ndiv / length(grouping) / n / nr)
}


monotoneEM <- function(x){
    mono <- 0
    for(i in 1:length(x)){
        if(any(x[[i]][-1] - x[[i]][-length(x[[i]])] < -.Machine$double.eps)){
            mono <- mono + 1
            print(i)
        }
    }
    mono
}

rObs <- function(x){
    temp <- x[-length(x)] == x[length(x)]
    if(sum(temp) == 0){
        tmax <- max(x[-length(x)][which(x[length(x)] - x[-length(x)] >0)])
        temp <- x[-length(x)] == tmax
    }
    if(sum(temp) > 1) temp[sample(which(temp) ,sum(temp)-1 )] <- FALSE
    which(temp)
}
