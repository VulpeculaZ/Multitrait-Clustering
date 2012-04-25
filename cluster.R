postRcpp <- cxxfunction(signature(data ="list",  P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/likelihood.cpp"), collapse = "\n" ))

cvRcpp <- cxxfunction(signature(data ="list", cvData="list",P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/cvlikelihood.cpp"), collapse = "\n" ))


##' <description>
##' Cluster observations with multiple traits
##' <details>
##' @title Multi-Traits Cluster
##' @param seqData Sequence Data.
##' @param P Permutation probability matrix.
##' @param bf Base Frequency.
##' @param obsHost The vector of hosts for the observations.
##' @param distMat Distance Matrix.
##' @param M Number of Clusters.
##' @param initP Initial Probability.
##' @param maxIter Maximum number of iterations.
##' @param tol Tolerance
##' @param method Which method to use? "rlog" or "cpp"
##' @param sparse Use sparse matrix for distance matrix clades?
##' @param model The probabilistic model for the hosts.
##' @return A list. $Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
mtCluster <- function(seqData, P, obsHost, distMat, initP,  bf=c(0.25,0.25,0.25,0.25), maxIter=200, tol=1e-4, method = "cpp", sparse = FALSE, model="geom", CV=NULL){
    M <- dim(initP)[2]
    n <- dim(initP)[1]
    logEzHost <- Ez <- initP
    logLike.old <- NULL
    p <- colSums(initP) / n
    ## for sequence
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    P <- as.matrix(P)
    ## for hosts
    hostLvl <- max(distMat)
    hostNum <- dim(distMat)[1]
    q <- matrix(, hostNum, M)
    hostClusters <- matrix(NA, hostNum, M)
    obsDist <- distMat[obsHost,]
    lambda <- rep(0.5, M)
    for(j in 1:maxIter){
        ## host:
        if(model == "geom"){
            for(i in 1:M){
                poisObs <- matrix(,n, hostNum)
                poisP <- dgeom(0:hostLvl, lambda[i], log=TRUE)
                for(k in 1:hostNum){
                    poisObs[,k] <- poisP[obsDist[,k]+1]
                }
                qMat <- sweep(poisObs, 1, Ez[,i], FUN ="*")
                q[,i] <- colSums(qMat)
                q[,i] <- exp(q[,i] - logsumexp(q[,i]))
                lambda[i] <- 1 / (1 + sum(exp(qMat) * obsDist) / sum(exp(qMat)))
                poisP <- dgeom(0:hostLvl, lambda[i], log=TRUE)
                for(k in 1:hostNum){
                    poisObs[,k] <- poisP[obsDist[,k]+1]
                }
                qMat <- sweep(poisObs, 2, log(q[,i]), FUN ="+")
                logEzHost[,i] <- apply(qMat, 1 , logsumexp)
            }
        }
        ## sequence:
        if(method=="rlog"){
            logec <- logPostR(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }
        if(method=="cpp"){
            logec <- postRcpp(seqData, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        }

        logEz <- sweep(logec + logEzHost, 2, log(p), FUN="+")
        logRowSumsEz <- apply(logEz, 1, logsumexp)
        logLike <- sum(logRowSumsEz)
        Ez <- exp(sweep(logEz, 1, logRowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            return(list(Ez = Ez, Ez.old = Ez.old, p = p, iter = j, convergence = convergence, logLike = logLike, lambda = lambda, q = q))
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
    if(!is.null(CV)){
        n <- length(CV$obsHost)
        logEzHost <- matrix(,n,M)
        cvDist <- distMat[CV$obsHost,]
        if(model == "geom"){
            for(i in 1:M){
                poisObs <- matrix(,n, hostNum)
                poisP <- dgeom(0:hostLvl, lambda[i], log=TRUE)
                for(k in 1:hostNum){
                    poisObs[,k] <- poisP[cvDist[,k]+1]
                }
                qMat <- sweep(poisObs, 2, log(q[,i]), FUN ="+")
                logEzHost[,i] <- apply(qMat, 1 , logsumexp)
            }
        }
        logec <- cvRcpp(seqData, cvData=CV$seq, P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        logEz <- sweep(logec + logEzHost, 2, log(p), FUN="+")
        logRowSumsEz <- apply(logEz, 1, logsumexp)
        logLikeCV <- sum(logRowSumsEz)
        return(list(Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike,lambda = lambda, q = q, logLikeCV=logLikeCV))
    }

    list(Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike,lambda = lambda, q = q)
}

##' <description>
##' A function for doing cross-validation to determine the number of clusters.
##' <details>
##' @title Cross-validation for clustering.
##' @param dataList Data is a list to be passed to clustering function.
##' @param n Number of observations.
##' @param maxCluster The maximum number of clusters. The cross-validation will try
##' from 1 to maxCluster number of clusters.
##' @param mCV The number of Monte Carlo drops for CV.
##' @param beta The ratio of learning samples.
##' @return A vector of mean log-likelihood for different number of clusters.
##' @author Ziqian Zhou
cvCluster <- function(dataList, n, maxCluster, mCV, beta){
    nCV <- ceiling(n * beta)
    logLikeCV <- matrix(0, mCV, length(maxCluster))
    for(i in 1:mCV){
        trSample <- sample(1:n, nCV)
        trSeq <- dataList$obsSeq[trSample]
        cvSeq <- dataList$obsSeq[-trSample]
        attributes(trSeq)$contrast <- attributes(cvSeq)$contrast <- attributes(dataList$obsSeq)$contrast
        attributes(trSeq)$nc <- attributes(cvSeq)$nc <- attributes(dataList$obsSeq)$nc
        cvData <- list(obsHost=dataList$obsHost[-trSample], seq=cvSeq)
        for(j in maxCluster){
            initP <- randomPMat(0.05, j, nCV)
            templl <- mtCluster(trSeq, dataList$P, dataList$obsHost[trSample], distMat = dataList$distMat , initP = initP, CV=cvData)
            logLikeCV[i,which(j == maxCluster)] <- templl$logLikeCV
        }
    }
    logLikeCV <- logLikeCV / mCV
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
##' @param obsHost A vector for the specis of hosts for the samples
##' @param distMat The distance matrix between hosts.
##' @param M Number of Clusters
##' @param initP Initial Probability
##' @param maxIter Maximum number of iterations
##' @param model. Available probabilistic models are geometric and poission(?).
##' @param tol Tolerance for convergence.
##' @return A list. Ez is the fitted value for clustering probability.
##' @author Ziqian Zhou
hostCluster <- function(obsHost, distMat, initP, maxIter=100, model="geom", tol=1e-4){
    M <- dim(initP)[2]
    logEzHost <- Ez <- initP
    logLike.old <- NULL
    n <- length(obsHost)
    p <- colSums(initP) / n
    hostLvl <- max(distMat)
    hostNum <- max(obsHost)
    q <- matrix(, hostNum, M)
    hostClusters <- matrix(NA, hostNum, M)
    obsDist <- distMat[obsHost,]
    lambda <- rep(0.5, M)
    poisP <- NA
    ## poisP <- matrix(dpois(x = rep(0:hostLvl, M) 0.5), M, hostLvl)
    for(j in 1:maxIter){
        ## Get the multinomial distribution for the hosts.
        ## get lambda for the transition probability
        if(model == "geom"){
            for(i in 1:M){
                ## M step:
                poisObs <- matrix(,n, hostNum)
                ## lambda: M * hostNum matrix, with estimated lambda for poisson
                ## lambda[i,] <- apply(wX, 2, sum) / p[i] / n
                ## E step
                ## poisP: matrix for lookup the pmf for cluster i.
                ## A hostLvl * hostNum matrix
                poisP <- dgeom(0:hostLvl, lambda[i], log=TRUE)
                for(k in 1:hostNum){
                    poisObs[,k] <- poisP[obsDist[,k]+1]
                }
                qMat <- sweep(poisObs, 1, Ez[,i], FUN ="*")
                q[,i] <- colSums(qMat)
                q[,i] <- exp(q[,i] - logsumexp(q[,i]))
                lambda[i] <- 1 / (1 + sum(exp(qMat) * obsDist) / sum(exp(qMat)))
                poisP <- dgeom(0:hostLvl, lambda[i], log=TRUE)
                for(k in 1:hostNum){
                    poisObs[,k] <- poisP[obsDist[,k]+1]
                }
                qMat <- sweep(poisObs, 2, log(q[,i]), FUN ="+")
                logEzHost[,i] <- apply(qMat, 1 , logsumexp)
            }
        }
        Ez <- sweep(logEzHost, 2, log(p), FUN="+")
        rowSumsEz <- apply(Ez, 1, logsumexp)
        logLike <- sum(rowSumsEz)
        Ez <- exp(sweep(Ez, 1, rowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        if(any(is.na(Ez))){
            warning("NA is produce in the last iteration step. Convergence failed.")
            return(list(Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, lambda = lambda, q = q))
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
         Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, lambda = lambda, q = q)
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


