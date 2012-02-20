##' <description>
##' Cluster observations with multiple traits
##' <details>
##' @title Multi-Traits Cluster
##' @param seqData Sequence Data
##' @param P Permutation probability matrix
##' @param bf Base Frequency
##' @param distMat Distance Matrix
##' @param M Number of Clusters
##' @param initP Initial Probability
##' @param maxIter Maximum number of iterations
##' @param tol Tolerance
##' @param method Which method to use? "rlog" or "cpp"
##' @param sparse Use sparse matrix for distance matrix clades?
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
mtCluster <- function(seqData, P, bf=c(0.25,0.25,0.25,0.25),  distMat, M, initP, maxIter, tol, method = "cpp", sparse = FALSE){
    distData <- designObsClade(distMat, sparse = FALSE)
    y <- distData$y
    X <- distData$X
    wMat <- combp(X, initP)
    Ez <- initP
    p <- rep(1/M, M)
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    P <- as.matrix(P)
    x.old <- 1

    for(j in 1:maxIter){
        ## M step:
        if(sparse == FALSE){
            mCluster <- cluster.nnls(y, X, w = wMat, sparse = sparse)
        }
        sigma.new <- apply(mCluster$x, 2, sd)
        ## E step:
        if(method=="rlog") ec <- logPostR(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        if(method=="cpp") ec <- postRcpp(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        for(i in 1:M){
            Ez[,i] <- p[i] * dnorm(mCluster$x[,i], sd=sigma.new[i]) * ec[,i]
        }
        Ez <- sweep(Ez, 1, rowSums(Ez), "/")
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        wMat <- combp(X, Ez)
        if(is.na(sum((x.old - mCluster$x)^2) < tol)) break
        if(sum((x.old - mCluster$x)^2) < tol) break
        x.old <- mCluster$x
    }
    list(fitted = mCluster$fitted, Ez = Ez, p = p, iter = j, convergence = sum((x.old - mCluster$x)^2), x=mCluster$x)
}
