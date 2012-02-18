##' <description>
##' Use EM algorithm to calculate the probability of for clustering.
##' <details>
##' @title
##' @param tree
##' @param seqData
##' @param P
##' @param bf Base frequency.
##' @param init Initial probability.
##' @param iter Number of iterations.
##' @param method Whether to use R or C implement, whether to use logarithm is probability calculation. (Temporary)
##' @return A matrix of probabilty.
##' @author Ziqian Zhou
seqEM <- function(tree, seqData, P, bf=c(0.25,0.25,0.25,0.25), init, iter=50, method="cpp"){
    nrs <- as.integer(length(seqData[[1]]))
    ncs <- as.integer(attr(seqData,"nc"))
    contrast <- attr(seqData, "contrast")
    ncos <- dim(contrast)[1]
    ec <- init
    P <- as.matrix(P)
    for(i in 1:iter){
        if(method=="r")  ec <- postR(seqData[tree$tip.label], P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=ec)
        if(method=="rlog") ec <- logPostR(seqData[tree$tip.label], P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=ec)
        if(method=="cpp") ec <- postRcpp(seqData[tree$tip.label], P=P, contrast=contrast, nrs=nrs, ncs=ncs, ncos=ncos, bf=bf, ecps=ec)
        ec <- ec/rowSums(ec)
    }
    ec
}

postRcpp <- cxxfunction(signature(data ="list", P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=srcLike)



postR <- function(seqData, P, contrast, nrs, ncs, ncos, bfs=c(0.25,0.25,0.25,0.25), ecps){
    ncluster=dim(ecps)[2]
    n <- length(seqData)
    post <- matrix(1, nrs, ncs)
    tmpPost <- post
    lookup <- contrast %*% P
    ecn <- matrix(nrow=n,ncol=ncluster)

    for(l in 1:ncluster){
        post <- matrix(1, nrs, ncs)
        for(i in 1:n){
            tmpSeq <- seqData[[i]]
            for(j in 1:nrs){
                tmpPost[j,] <- lookup[tmpSeq[j],]
            }
            post = post * (tmpPost ^ ecps[i,l])
            ecn[i,l] <- mean(post)
        }

        for(i in 1:n){
            tmpSeq <- seqData[[i]]
            for(j in 1:nrs){
                tmpPost[j,] <- lookup[tmpSeq[j],]
            }
            posti <- post /  tmpPost^ecps[i,l]
            posti <- posti / rowSums(posti)
            ecn[i,l] <- prod((posti * tmpPost) %*% bfs)
        }
    }
    return(ecn)
}


logPostR <- function(seqData, P, contrast, nrs, ncs, ncos, bfs=c(0.25,0.25,0.25,0.25), ecps){
    ncluster=dim(ecps)[2]
    n <- length(seqData)
    post <- matrix(1, nrs, ncs)
    tmpPost <- log(post)
    lookup <- log(contrast %*% P)
    ecn <- matrix(nrow=n,ncol=ncluster)
    for(l in 1:ncluster){
        post <- matrix(0, nrs, ncs)
        for(i in 1:n){
            tmpSeq <- seqData[[i]]
            for(j in 1:nrs){
                tmpPost[j,] <- lookup[tmpSeq[j],]
            }
            post = post + (tmpPost * ecps[i,l])
            ecn[i,l] <- mean(post)
        }
        for(i in 1:n){
            tmpSeq <- seqData[[i]]
            for(j in 1:nrs){
                tmpPost[j,] <- lookup[tmpSeq[j],]
            }
            posti <- exp(post - (tmpPost*ecps[i,l]))
            posti <- posti / rowSums(posti)
            ecn[i,l] <- prod((posti * exp(tmpPost)) %*% bfs)
        }
    }
    return(ecn)
}


simseqEM <- function(tc, iter=100, seqLength=100, initRange = 0.05, edgeLength, P, method="cpp"){
    nTips <- length(tc$tip.label)
    if(is.null(tc$edge.length)){
        inner <- which(tc$edge[,2]>nTips)
        tc$edge.length <- rep(0.1, length(tc$edge[,1]))
        tc$edge.length[inner] <- edgeLength
    }
    result <- list()
    for(i in 1:iter){
        init <- runif((tc$Nnode-1)*nTips, -initRange, initRange)
        ecps <- matrix(1 / (tc$Nnode-1) +init , nrow = nTips, ncol = tc$Nnode-1)
        ecps <- ecps / rowSums(ecps)
        sc <- simSeq(tc, l=seqLength)
        result[[i]]  <- seqEM(tree=tc, seqData=sc, P=P, init=ecps, iter=50, method=method)
    }
    result
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



##' <description>
##'
##' <details>
##' @title
##' @return
##' @author Ziqian Zhou
##' @param nodesL
##' @param seqList
seqPost <- function(nodesL, seqList)
{
    result <- list()
    nTips <- length(seqList)
    nr <- attr(seqList, "nr")
    for(i in 2:length(nodesL)){
        result[[i-1]] <- rep(1, nTips)
        for(j in 1:nTips){
            for(k in 1:nr){
                result[[i-1]][j] <- result[[i-1]][j] * nodesL[[i]][k,seqList[[j]][k]]
            }
        }
    }
    result
}


pcl <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1,
    rate = 1, model=NULL, ...)
{
    call <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix")
    existing <- match(pmla, names(extras))
    wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], parent.frame()))
    llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], parent.frame()))

    if(class(tree)!="phylo") stop("tree must be of class phylo")
    if (is.null(attr(tree, "order")) || attr(tree, "order") ==
        "cladewise")
        tree <- phangorn:::reorderPruning(tree)
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(is.null(tree$edge.length)) stop("tree must have edge weights")
    if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    nr <- attr(data, "nr")
    ## length of sequence?
    type <- attr(data,"type")
    if(type=="AA" & !is.null(model)){
        model <- match.arg(model, c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
### Not dealing with AA right now.
### Need to be done later.
### getModelAA(model, bf=is.null(bf), Q=is.null(Q))
    }

    if (is.null(bf))
        bf <- rep(1/length(levels), length(levels))
    if (is.null(Q))
        Q <- rep(1, length(levels) * (length(levels) - 1)/2)
    m <- 1
    eig <- phangorn:::edQt(bf = bf, Q = Q)
    w <- rep(1/k, k)
    if (inv > 0)
        w <- (1 - inv) * w
    if (wMix > 0)
        w <- wMix * w
    g <- phangorn:::discrete.gamma(shape, k)
    if (inv > 0)
        g <- g/(1 - inv)
    g <- rate * g
## #    INV <- lli(data, tree)
##     INV <- Matrix(lli(data, tree), sparse=TRUE)
## #    ll.0 <- INV %*% (bf * inv)
##     ll.0 <- as.matrix(INV %*% (bf * inv))
##     if(wMix>0) ll.0 <- ll.0 + llMix
    resll <- matrix(0, nr, k)
    while (m <= k) {
        resll <- llan(data, tree, bf = bf, g = g[m], Q = Q, eig = eig, assign.dat = FALSE, ...)
        m = m + 1
    }
    resll
}



##' <description>
##'
##' <details>
##' This function is adapted from ll function from phangorn package.
##' @title
##' @param dat1
##' @param tree
##' @param bf
##' @param g
##' @param Q
##' @param eig
##' @param assign.dat
##' @param ...
##' @return A list. The likelihood of all nodes.
##' @author Ziqian Zhou
llan <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1,
    Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE, ...)
{
#    if(is.null(attr(tree,"order")) || attr(tree,"order")=="cladewise")tree <- reorderPruning(tree)
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1
#    dat[1:q] = dat1[tree$tip.label]
    if (is.null(eig)) eig = edQt(bf = bf, Q = Q)
    el <- tree$edge.length
    P <- phangorn:::getP(el, eig, g)
    nr <- as.integer(attr(dat1,"nr"))
    nc <- as.integer(attr(dat1,"nc"))
    node = as.integer(node-min(node))
    edge = as.integer(edge-1L)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1L)
    contrast = attr(dat1, "contrast")
    nco = as.integer(dim(contrast)[1])
    res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
    result <- res
    return(list(result=result, P=P, contrast=contrast))
}
