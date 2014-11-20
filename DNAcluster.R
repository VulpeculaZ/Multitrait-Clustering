library(phangorn)
library(inline)
library(RcppArmadillo)
library(cluster)

## library(combinat)
## library(parallel)
## library(nnet)

emSeqRcpp <- cxxfunction(signature(data ="list", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", nts="integer", ecps="matrix", weightS="integer",  MachineEps="double"), plugin="RcppArmadillo", body=paste( readLines("./Source/emSeq.cpp"), collapse = "\n" ))



.nIdenVec <- function(x, y, weight){
    sum((unlist(x) !=y) * weight)
}

## First M diagnal are for theta, the last M-1 are for p.
varTheta <- function(cluster, obs){
    theta <- cluster$theta
    y <- cluster$y
    p <- cluster$p
    Ez <- cluster$Ez
    nY <- length(attr(obs, "index"))
    weightS <- attr(obs, "weight")
    logTheta <- log(theta)
    log13Theta <- log(1 - 3 * theta)
    M <- length(p)
    ## obs Fisher information
    obsDiff <- matrix(, dim(Ez)[1], dim(Ez)[2])
    for(i in 1:M){
        obsDiff[,i] <- mapply(.nIdenVec, x=obs, MoreArgs = list(y = y[,i], weight = weightS))
    }
    ## Ix|y
    plogTheta <- sweep(obsDiff, 2, theta, FUN="/") - sweep(3*(nY-obsDiff), 2, (1-3*theta), FUN="/")
    varZ <- Ez * (1-Ez)
    varStheta <- colSums(plogTheta ^ 2 * varZ)
    varZp2 <- sweep(varZ, 2, p^2 , "/")
    plogObs <- Ez * (sweep(obsDiff, 2, theta, FUN="/") - sweep(3*(nY-obsDiff), 2, (1-3*theta), FUN="/"))
    if(M>2){
        covZim <- sweep(Ez[,1:(M-1)], 1, -Ez[,M], "*")
        varSp <- colSums(sweep(varZp2[,1:(M-1)], 1, varZp2[,M], "+") + 2 * sweep(covZim, 2, p[1:(M-1)], FUN = "/") / p[M])
    }
    else{
        covZim <- matrix(-Ez[,1] * Ez[,2], nrow = dim(Ez)[1])
        varSp <- sum(varZp2[,1] + varZp2[,2] + 2 * covZim / p[1] /p[M])
    }
    Ixy <- matrix(, M + M -1, M + M -1)
    diag(Ixy) <- c(varStheta,varSp)
    ## Cov of S_{c,theta_i} and S_{c,theta_j}
    for(i in 2 : M){
        for(j in 1:(i-1)){
            Ixy[i, j] <- sum(plogTheta[,i] * plogTheta[,j] * Ez[,i] * Ez[,j])
        }
    }
    ## Cov of S_{c,pi} and S_{c,pj}
    if(M > 2){
        for(i in 2: (M - 1)){
            for(j in 1:(i-1)){
                Ixy[i + M, j+M] <- sum(- Ez[,i] * Ez[,j] / (p[i]* p[j]) - covZim[,i] / (p[i] * p[M]) - covZim[,j] / (p[j] * p[M]) + varZ[,M] / (p[M]^2))
            }
        }
    }
    ## Cov of S_{c,pi} and S_{c,theta_j}
    for(i in 1:(M-1)){
        Ixy[i+M, M] <- sum(plogTheta[,M] * (covZim[,i] / p[i] - varZ[,M] / p[M] ))
        for(j in 1:(M-1)){
            if(j == i){
                Ixy[i + M,j] <- sum(plogTheta[,j] * (varZ[,j] / p[i] - covZim[,i] / p[M]))
            }
            else{
                Ixy[i + M,j] <- sum(plogTheta[,j] * (- Ez[,j] * Ez[,i] / p[i] - covZim[,i] / p[M]))
            }
        }
    }
    Ixy[upper.tri(Ixy)] <- t(Ixy)[upper.tri(Ixy)]
    ## Ic
    Ezp <- sweep(Ez, 2, p^2, FUN = "/")
    if(M>2){
        Ezpim <- sweep(Ezp[,1:(M-1)], 1, Ezp[,M], "+")
    }
    else{
        Ezpim <- matrix(Ezp[,1] + Ezp[,M], nrow = dim(Ez)[1])
    }
    Ic <- diag(c(colSums(Ez * ( sweep(obsDiff, 2, theta^2, FUN="/")  + sweep(9*(nY-obsDiff), 2, (1-3*theta)^2, FUN="/"))), colSums(Ezpim)))
    varScore <- Ic - Ixy
    solve(varScore)
}


fNbi <- function(seq, ezil, contrast){
    tmpSeq <- contrast[seq,]
    tmpSeq * ezil
}

fLogp <- function(seq, yl, nrs, log1N3theta, logtheta){
    nSame <- sum(yl == seq)
    nSame* log1N3theta + (nrs - nSame) * logtheta
}


EMinit <- function(distMat, k, diffP=0.05){
    n <- dim(distMat)[1]
    clusterLab <- pam(x = distMat, k=k, diss = FALSE, cluster.only=TRUE)
    initP <- matrix(0, nrow = n, ncol = k)
    initP[cbind(1:n, clusterLab)] <- diffP
    initP <- initP + (1-diffP) / k
    initP
}

sumLogLike <- function(x, M, logLike){
    sum(logLike[cbind(1:M, x)])
}

logsumexp <- function(x){
   xmax <- which.max(x)
   log1p(sum(exp(x[-xmax]-x[xmax]), na.rm = TRUE)) + x[xmax]
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
DNACluster <- function(seqData, bf=c(0.25,0.25,0.25,0.25), initP, maxIter=200, tol=1e-4, method = "cpp", CV=NULL, P=NULL, nSame=0){
    M <- dim(initP)[2]
    n <- dim(initP)[1]
    Ez <- initP
    p <- colSums(initP) / n
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    nts <- length(attr(seqData, "index"))
    weightS <- attr(seqData, "weight")
    ##P <- list()
    logLike.old <- NULL
    distDNA <- matrix(,n,n)
    y <- NULL
    theta <- NULL
    for(j in 1:maxIter){
        ## E step:
        if(method=="EMseq"){
            emList <- EMseqR(seqData, contrast=contrast, nrs=nrs, ncs=ncs, ecps=Ez)
            logec <- emList$ecps
            y <- emList$y
            theta <- emList$theta
        }
        if(method =="cpp"){
            browser(expr = !is.matrix(Ez))
            ## save(seqData, contrast, nrs, ncs, ncos, Ez, file = "failData.RData")
            emList <- emSeqRcpp(data=seqData, contrast=contrast, nrs=nrs, ncs =ncs, ncos=ncos, nts =nts, ecps=Ez, weightS = weightS,  MachineEps=.Machine$double.eps)
            logec <- emList$ecps
            y <- emList$y
            theta <- emList$theta
        }
        logEz <- sweep(logec, 2, log(p), FUN="+")
        logRowSumsEz <- apply(logEz, 1, logsumexp)
        logLike <- sum(logRowSumsEz)
        Ez <- exp(sweep(logEz, 1, logRowSumsEz, "-"))
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        ## if(any(is.na(Ez))){
        ##     warning("NA is produce in the last iteration step. Convergence failed.")
        ##     return(list(Ez = Ez, Ez.old = Ez.old, p = p, iter = j, convergence = convergence, logLike = logLike, P=P, y=y))
        ## }
        if(any(is.na(Ez))){
            Ez[which(is.na(Ez))] <- 0
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
                if(logLike.init.increase <0) logLike.init.increase <- logLike - logLike.old
                if(logLike - logLike.old < logLike.init.increase * tol && logLike - logLike.old >= 0) break
                if(abs(logLike - logLike.old) < .Machine$double.eps ) break
            }
        }
        logLike.old <- logLike
    }
    if(!is.null(CV)){
        n <- length(CV$seq)
        logEzHost <- matrix(,n,M)
        ## logec <- cvRcpp(seqData, cvData=CV$seq, PList=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=Ez)
        logec <- seqCV(cvData=CV$seq,  contrast=contrast, nrs=nrs, ncs=ncs, ecps=Ez, y =y, theta= theta)
        logEz <- sweep(logec, 2, log(p), FUN="+")
        logRowSumsEz <- apply(logEz, 1, logsumexp)
        logLikeCV <- sum(logRowSumsEz)
        return(list(Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, logLikeCV=logLikeCV))
    }
    list(Ez = Ez, p = p, iter = j, convergence = convergence, logLike = logLike, P=P, y=y, theta = theta)
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
cvCluster <- function(dataList, n, maxCluster, mCV, beta, initPList, method="DNA"){
    nCV <- ceiling(n * beta)
    logLikeCV <- matrix(NA, mCV, maxCluster)
    for(i in 1:mCV){
        trSample <- sample(1:n, nCV)
        trSeq <- dataList$obsSeq[trSample]
        cvSeq <- dataList$obsSeq[-trSample]
        attributes(trSeq)$contrast <- attributes(cvSeq)$contrast <- attributes(dataList$obsSeq)$contrast
        attributes(trSeq)$nc <- attributes(cvSeq)$nc <- attributes(dataList$obsSeq)$nc
        attributes(trSeq)$weight <- attributes(cvSeq)$weight <- attributes(dataList$obsSeq)$weight
        attributes(trSeq)$index <- attributes(cvSeq)$index <- attributes(dataList$obsSeq)$index


        if(method =="mt"){
            cvData <- list(obsHost=dataList$obsHost[-trSample], seq=cvSeq)
            for(j in 1:maxCluster){
                initP <- initPList[[j]][trSample,]
                templl <- mtCluster(trSeq, dataList$obsHost[trSample], distMat = dataList$distMat , initP = initP, CV=cvData)
                try(logLikeCV[i,j] <- templl$logLikeCV)
            }
        }
        if(method == "DNA"){
            cvData <- list(seq = cvSeq)
            for(j in 1:maxCluster){
                initP <- initPList[[j]][trSample,]
                if(j == 1) initP <- matrix(1, length(trSample),1)
                templl <- DNACluster(trSeq, initP = initP, CV=cvData)
                try(logLikeCV[i,j] <- templl$logLikeCV)
            }
        }
    }
    logLikeCV
}




## pcl <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1,
##     rate = 1, model=NULL, ...)
## {
##     call <- match.call()
##     extras <- match.call(expand.dots = FALSE)$...
##     pmla <- c("wMix", "llMix")
##     existing <- match(pmla, names(extras))
##     wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], parent.frame()))
##     llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], parent.frame()))

##     if(class(tree)!="phylo") stop("tree must be of class phylo")
##     if (is.null(attr(tree, "order")) || attr(tree, "order") ==
##         "cladewise")
##         tree <- phangorn:::reorderPruning(tree)
##     if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
##     if(is.null(tree$edge.length)) stop("tree must have edge weights")
##     if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
##     levels <- attr(data, "levels")
##     weight <- attr(data, "weight")
##     nr <- attr(data, "nr")
##     ## length of sequence?
##     type <- attr(data,"type")
##     if(type=="AA" & !is.null(model)){
##         model <- match.arg(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
## ### Not dealing with AA right now.
## ### Need to be done later.
## ### getModelAA(model, bf=is.null(bf), Q=is.null(Q))
##     }

##     if (is.null(bf))
##         bf <- rep(1/length(levels), length(levels))
##     if (is.null(Q))
##         Q <- rep(1, length(levels) * (length(levels) - 1)/2)
##     m <- 1
##     eig <- phangorn:::edQt(bf = bf, Q = Q)
##     w <- rep(1/k, k)
##     if (inv > 0)
##         w <- (1 - inv) * w
##     if (wMix > 0)
##         w <- wMix * w
##     g <- phangorn:::discrete.gamma(shape, k)
##     if (inv > 0)
##         g <- g/(1 - inv)
##     g <- rate * g
## ## #    INV <- lli(data, tree)
## ##     INV <- Matrix(lli(data, tree), sparse=TRUE)
## ## #    ll.0 <- INV %*% (bf * inv)
## ##     ll.0 <- as.matrix(INV %*% (bf * inv))
## ##     if(wMix>0) ll.0 <- ll.0 + llMix
##     resll <- matrix(0, nr, k)
##     while (m <= k) {
##         resll <- llan(data, tree, bf = bf, g = g[m], Q = Q, eig = eig, assign.dat = FALSE, ...)
##         m = m + 1
##     }
##     resll
## }



## ##' <description>
## ##'
## ##' <details>
## ##' This function is adapted from ll function from phangorn package.
## ##' @title
## ##' @param dat1
## ##' @param tree
## ##' @param bf
## ##' @param g
## ##' @param Q
## ##' @param eig
## ##' @param assign.dat
## ##' @param ...
## ##' @return A list. The likelihood of all nodes.
## ##' @author Ziqian Zhou
## llan <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1,
##     Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE, ...)
## {
## #    if(is.null(attr(tree,"order")) || attr(tree,"order")=="cladewise")tree <- reorderPruning(tree)
##     q = length(tree$tip.label)
##     node <- tree$edge[, 1]
##     edge <- tree$edge[, 2]
##     m = length(edge) + 1
## #    dat[1:q] = dat1[tree$tip.label]
##     if (is.null(eig)) eig = edQt(bf = bf, Q = Q)
##     el <- tree$edge.length
##     P <- phangorn:::getP(el, eig, g)
##     nr <- as.integer(attr(dat1,"nr"))
##     nc <- as.integer(attr(dat1,"nc"))
##     node = as.integer(node-min(node))
##     edge = as.integer(edge-1L)
##     nTips = as.integer(length(tree$tip))
##     mNodes = as.integer(max(node) + 1L)
##     contrast = attr(dat1, "contrast")
##     nco = as.integer(dim(contrast)[1])
##     res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
##     result <- res
##     return(list(result=result, P=P, contrast=contrast))
## }
