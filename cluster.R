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


##' <description>
##' Cluster observations with distance matrix
##' <details>
##' @title Multi-Traits Cluster
##' @param distMat Distance Matrix
##' @param M Number of Clusters
##' @param initP Initial Probability
##' @param maxIter Maximum number of iterations
##' @param tol Tolerance
##' @param sparse Use sparse matrix for distance matrix clades?
##' @param tol Tolerance,
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
DistCluster <- function(P, bf=c(0.25,0.25,0.25,0.25),  distMat, M, initP, maxIter, tol=1e-4, method = "cpp", sparse = FALSE){
    distData <- designObsClade(distMat, sparse = FALSE)
    y <- distData$y
    X <- distData$X
    wMat <- combp(X, initP)
    Ez <- initP
    p <- rep(1/M, M)
    x.old <- 1
    logLike.old <- NULL
    n <- dim(distMat)[1]
    combn2 <- t(combn(n,2))
    for(j in 1:maxIter){
        ## M step:
        if(sparse == FALSE){
            mCluster <- cluster.nnls(y, X, w = wMat, sparse = sparse)
        }
        ## ME
        ## sigma.new <- sqrt(colSums(mCluster$x^2 * Ez) / colSums(Ez))
        ## LS
        sigma.new <- sqrt(mCluster$sigma / colSums(wMat))
        logf2 <- apply(rbind(mCluster$res, sigma.new), 2, function(x) dnorm(x[1:(length(x)-1)], sd = x[length(x)], log=TRUE))
        Ez.new <- matrix(0, n, M)
        matlogf <- list()
        for(i in 1:M){
            matlogf[[i]] <- matrix(0, n, n)
            matlogf[[i]][combn2] <- logf2[,i]
            matlogf[[i]] <- matlogf[[i]] + t(matlogf[[i]])
        }
        if(any(is.na(matlogf[[1]]))){
            browser()
        }
        Ez <-  mcEz(matlogf, p = p)
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
        if(j == 1){
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
##' Cluster observations with host information
##' <details>
##' @title Cluster Using Host information
##' @param host A vector for the specis of hosts for the samples
##' @param P A probability for changing the host
##' @param hostMat Matrix of the host tree structure
##' @param M Number of Clusters
##' @param initP Initial Probability
##' @param maxIter Maximum number of iterations
##' @param tol Tolerance
##' @return A list. Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
hostCluster <- function(obsHost, distMat, hostMat, M, initP, maxIter, method="poisson", tol=1e-4){
    Ez <- initP
    p <- colSums(initP)
    logLike.old <- NULL
    n <- length(obsHost)
    hostLvl <- max(hostMat[,1])
    hostNum <- dim(hostMat)[2]
    qMat <- matrix(nrow = hostNum, ncol = M)
    hostClusters <- matrix(NA, hostNum, M)
    obsDist <- distMat[obsHost,]
    lambda <- matrix(, M, hostNum)
    for(j in 1:maxIter){
        ## M step:
        ## Get the multinomial distribution for the hosts.
                ## get lambda for the transition probability
        if(method == "poisson"){
            for(i in 1:M){
                poisObs <- matrix(,n, hostNum)
                wX <- sweep(obsDist, 1, Ez[,i], FUN="*")
                ## lambda: M * hostNum matrix, with estimated lambda for poisson
                lambda[i,] <- apply(wX, 2, sum) / p[i]
                ## poisP: matrix for lookup the pmf for cluster i.
                ## A hostLvl * hostNum matrix
                poisP <- vapply(lambda[i,], FUN=function(x) dpois(1:hostLvl, x), FUN.VALUE = 1:hostLvl)
                for(k in 1:hostNum){
                    poisObs <- poisP[obsDist[,k],k]
                }
                qMat[,i] <- colSum(sweep(log(poisObs), 1, Ez[,i]), FUN ="*")
                qMat[,i] <- exp(qMat[,i] - logsumexp(qMat[,i]))
                Ez[,i] <- rowSuma(sweep(poisObs,2,qMat[,i], FUN="*"))
            }
        }
        for(i in 1:hostNum){
            hostClusters[i,] <- colSums(sweep(Ez, 1,  (hostMat[,1] == i), FUN = "*"))
        }
        hostClusters <- sweep(hostClusters, 2, colSums(hostClusters), FUN ="/")

        ## E step:
        Ez.old <- Ez
        for(i in 1:M){
            probi <- dpois(1:hostLvl,lambda = lambda[i])
            obsP <- matrix(NA, n, hostNum)
            for(k in 1:hostLvl){
                obsP[obsDist == k] <- probi[k]
            }
            Ez[,i] <- obsP %*% t(hostClusters[,i]) * p[i]
        }
        logLike <- sum(log(rowSums(Ez)))
        Ez <- sweep(Ez, 1, rowSums(Ez), "/")
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            return(list(Ez = Ez, Ez.old = Ez.old, p = p, iter = j, convergence = convergence, logLike = logLike))
        }
        if(j == 1){
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
   log1p(sum(exp(x[-xmax]-x[xmax]), na.rm = TRUE)) + x[xmax]
}


mcEz <- function(matlogf, p, MC = "MC", iter = NULL){
    N <- dim(matlogf[[1]])[1]
    M <- length(p)
    logp <- log(p)
    Ez <- matrix(nrow = N, ncol=M)
    logLikeAll <- NA
    if(MC == "MC"){
        if(is.null(iter)) iter =  (M * N)^2 * floor(log(N * M))
        allIndex <- sample.int(n = M, size = iter * N, replace =TRUE)
        for(i in 1:iter){
            index <- allIndex[(1+ (i-1) * N): (i * N)]
            logLike <- 0
            for(j in 1:M){
                indexj <- (index == j)
                indexMat <- matrix(as.logical(outer(indexj,indexj)),N ,N)
                indexMat[!lower.tri(indexMat)] <- FALSE
                logLike <- logLike + sum(matlogf[[j]][indexMat]) + sum(indexj) * logp[j]
            }
            for(j in 1:N){
                temp <- logsumexp(c(Ez[j, index[j]], logLike))
                Ez[j, index[j]] <- temp
            }
        }
    }
    Ez
}


