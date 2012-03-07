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
##' @param tol Tolerance,
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
mtCluster <- function(seqData, P, bf=c(0.25,0.25,0.25,0.25),  distMat, M, initP, maxIter, tol=1e-4, method = "cpp", sparse = FALSE){
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
    logLike.old <- NULL
    for(j in 1:maxIter){
        ## M step:
        if(sparse == FALSE){
            mCluster <- cluster.nnls(y, X, w = wMat, sparse = sparse)
        }
        ## ME
        sigma.new <- sqrt(colSums(mCluster$x^2 * Ez) / colSums(Ez))
        ## LSE
        ## sigma.new <- sqrt(mCluster$sigma / colSums(wMat))
        ## E step:
        if(method=="rlog"){
            logec <- logPostR(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        if(method=="cpp"){
            logec <- postRcpp(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        ## ME
        for(i in 1:M){
            ## use dexp likelihood
            logfi <- (-mCluster$x[,i]) / sigma.new[i] - log(sigma.new[i])
            ## use normal likelihood
            ## logfi <- -(mCluster$x[,i]^2)
            Ez[,i] <- log(p[i]) + logfi + logec[,i]
        }
        ## LS
        ## logf2 <- apply(rbind(mCluster$res, sigma.new), 2, function(x) dnorm(x[1:(length(x)-1)], sd = x[length(x)], log=TRUE))
        ## for(i in 1:M){
        ##     logfi <- t(logf2[,i]) %*% X
        ##     Ez[,i] <- log(p[i]) + logfi + logec[,i]
        ## }
        logRowSumsEz <- apply(Ez, 1, logsumexp)
        logLike <- sum(logRowSumsEz)
        Ez <- exp(sweep(Ez, 1, logRowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        wMat <- combp(X, Ez)
        x.old <- mCluster$x
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            return(list(Ez = Ez, Ez.old = Ez.old, p = p, iter = j, convergence = convergence, logLike = logLike, x=mCluster$x))
        }
        if(j == 1) {
            convergence <- logLike
            Ez.old <- Ez
        }
        if(!is.null(logLike.old)){
            convergence <- c(convergence, logLike)
            if(j == 2){
                logLike.init.increase <- logLike - logLike.old
                Ez.old <- Ez
            }
            if(j > 2){
                Ez.old <- Ez
                if(logLike - logLike.old < logLike.init.increase * tol && logLike - logLike.old >= 0) break
            }
        }
        logLike.old <- logLike
    }
    list(#fitted = mCluster$fitted,
         Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, x=mCluster$x)
}


##' <description>
##' Cluster observations with DNA sequences
##' <details>
##' @title DNA cluster
##' @param seqData Sequence Data
##' @param P Permutation probability matrix
##' @param bf Base Frequency
##' @param M Number of Clusters
##' @param initP Initial Probability
##' @param maxIter Maximum number of iterations
##' @param tol Tolerance
##' @param method Which method to use? "rlog" or "cpp"
##' @param tol Tolerance,
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
DNACluster <- function(seqData, P, bf=c(0.25,0.25,0.25,0.25), M, initP, maxIter, tol=1e-4, method = "cpp"){
    Ez <- initP
    p <- rep(1/M, M)
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    P <- as.matrix(P)
    logLike.old <- NULL
    for(j in 1:maxIter){
        ## E step:
        if(method=="rlog"){
            logec <- logPostR(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        if(method=="cpp"){
            logec <- postRcpp(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        for(i in 1:M){
            Ez[,i] <- log(p[i]) + logec[,i]
        }
        logRowSumsEz <- apply(Ez, 1, logsumexp)
        logLike <- sum(logRowSumsEz)
        Ez <- exp(sweep(Ez, 1, logRowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            return(list(Ez = Ez, Ez.old = Ez.old, p = p, iter = j, convergence = convergence, logLike = logLike))
        }
        if(j == 1) {
            convergence <- logLike
            Ez.old <- Ez
        }
        if(!is.null(logLike.old)){
            convergence <- c(convergence, logLike)
            if(j == 2){
                logLike.init.increase <- logLike - logLike.old
                Ez.old <- Ez
            }
            if(j > 2){
                Ez.old <- Ez
                if(logLike - logLike.old < logLike.init.increase * tol && logLike - logLike.old >= 0) break
            }
        }
        logLike.old <- logLike
    }
    list(#fitted = mCluster$fitted,
         Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike)
}


logsumexp <- function(x){
   xmax <- which.max(x)
   log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
}
