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
##' @param tol Tolerance,
##' @param method Which method to use? "rlog" or "cpp"
##' @param sparse Use sparse matrix for distance matrix clades?
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
mtCluster <- function(seqData, P, bf=c(0.25,0.25,0.25,0.25),  distMat, M, initP, maxIter=200,  method = "cpp", sparse = FALSE, tol=1e-4){
    distData <- designObsClade(distMat, sparse = FALSE)
    y <- distData$y
    X <- distData$X
    n <- dim(X)[2]
    wMat <- combp(X, initP)
    Ez <- initP
    p <- rep(1/M, M)
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    P <- as.matrix(P)
    logLike.old <- NULL

    for(j in 1:maxIter){
        ## M step:
        if(sparse == FALSE){
            mCluster <- cluster.nnls(y, X, w = wMat, sparse = sparse)
            sigma.new <- sqrt(mCluster$sig)
            logf2 <- apply(rbind(mCluster$res, mCluster$sig), 2, function(x) dnorm(x[1:(length(x)-1)], sd = x[length(x)], log=TRUE))
            logf <- matrix(, n, M)
            for(k in 1:M){
                w4f <- t(mCluster$x[,k] * t(X))
                w4f <- sweep(w4f, 1, rowSums(w4f), "/")
                logf[,k] <- t(logf2[,k]) %*% w4f
            }
        }
        ## E step:
        if(method=="rlog"){
            logec <- logPostR(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        if(method=="cpp"){
            logec <- postRcpp(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        for(i in 1:M){
            Ez[,i] <- log(p[i]) + logf[,i] + logec[,i]
        }
        logRowSumsEz <- apply(Ez, 1, logsumexp)
        logLike <- sum(logRowSumsEz)
        Ez <- exp(sweep(Ez, 1, logRowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        wMat <- combp(X, Ez)
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            break
        }
        if(j == 1) convergence <- logLike
        if(!is.null(logLike.old)){
            convergence <- c(convergence, logLike)
            if(j == 2){
                logLike.init.increase <- logLike - logLike.old
            }
            if(j > 2){
                if(logLike - logLike.old < logLike.init.increase * tol && logLike - logLike.old > 0) break
            }
        }
        logLike.old <- logLike
    }
    list(#fitted = mCluster$fitted,
         Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, x=mCluster$x)
}
