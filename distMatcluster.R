cluster.em <- function(y, X, pMat, sparse = FALSE, maxIter = 100, tol=1E-12){
    M <- dim(pMat)[2]
    wMat <- combp(X, pMat)
    difference <- max(1, 100 * tol)
    Ez <- pMat
    p <- rep(1/M, M)
    for(j in 1:maxIter){
        ## M step:
        if(sparse == FALSE){
            mCluster <- cluster.nnls(y, X, w = wMat, sparse = sparse)
        }
        sigma.new <- sd(mCluster$x)
        ## E step:
        for(i in 1:M){
            Ez[,i] <- p[i] * dnorm(mCluster$x[,i], sd=sigma.new)
        }
        Ez <- sweep(Ez, 1, rowSums(Ez), "/")
        colSums.Ez <- colSums(Ez)
        p <- colSums.Ez / sum(colSums.Ez)
        wMat <- combp(X, Ez)
        if(j > 1){
            if(sum((x.old - mCluster$x)^2) < tol) break
        }
        x.old <- mCluster$x
    }
    list(fitted = mCluster$fitted, Ez = Ez, p = p, iter = j, convergence = sum((x.old - mCluster$x)^2))
}


##' ##' <description>
##' Fitting the cluster using non-negative least square.
##' <details>
##' @title cluster.nnls
##' @param y observed value
##' @param X Design matrix
##' @param pMat
##' @param sparse Is it a sparse matrix?
##' @return A list with residuals matrix res,  fitted value matrix fitted.
##' @author Ziqian Zhou
cluster.nnls <- function(y, X, w, sparse = FALSE){
    if(any(is.na(y)))
    {
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
    res <- fitted <- matrix(,dim(w)[1], dim(w)[2])
    x <- matrix(, dim(X)[2], dim(w)[2])
    if(sparse == FALSE){
        for(i in 1:dim(w)[2]){
            p <- w[,i]
            fit <- clade.nnls(y, X, p, sparse = sparse)
            res[,i] <- fit$residuals
            fitted[,i] <- fit$fitted
            x[,i] <- fit$x
        }
    }
    else{
    }
    list(res = res, fitted = fitted, x = x)
}

##' ##' <description>
##' Fitting the cluster using non-negative least square.
##' <details>
##' @title cluster.nnls
##' @param y observed value
##' @param X Design matrix
##' @param w Weight
##' @param sparse Is it a sparse matrix?
##' @return Fitted value
##' @author Ziqian Zhou
clade.nnls <- function(y, X, w=NULL, sparse = FALSE){
    ## need to use nnls() in R base
    ## or nnls for sparse matrices
    ## na.action
    ## LS solution
    if(sparse == FALSE){
        if(is.null(w)){
            fit <- nnls(X, y)
        }
        else{
            wsqrt <- sqrt(w)
            fit <- nnls(wsqrt*X, wsqrt*y)
            ## betahat = fit$coefficients
            ## RSS = sum(fit$residuals^2)
            ## print(paste("RSS:", RSS))
        }
        ## return(betahat)
    }
    else{
        ## X <- as(X,"sparseMatrix")
        ## fit <- MatrixModels:::lm.fit.sparse(y=y, x=X, w)
    }
    fit
}

combp <- function(X, pMat){
    n <- dim(pMat)[1]
    wMat <- matrix(, n*(n-1) / 2, dim(pMat)[2])
    for(i in 1:dim(pMat)[2]){
        wMat[,i] <- apply(X, 1, function(x) prod(pMat[as.logical(x),i]))
    }
    wMat
}

##' <description>
##' Create the design matrix and format observed values for least square criteria.
##' The dimension of the design matrix will be: n(n-1)(m^2)/2 by nm+m.
##' Observed
##' <details>
##' @title designObsCluster
##' @param n number of observation
##' @param yMat Input observed distances. Need to be in the form of
##' distance matrix.
##' @param sparse Will the returned matrix be sparse? Need to be a logical value.
##' @return A design Matrix representing the cluster by a tree structure. And a vector of observed distance.
##' @author Ziqian Zhou
designObsClade <- function( yMat, sparse=FALSE){
    if(!isSymmetric(yMat)) stop("the distance matrix must be symmetric")
    n <- dim(yMat)[1]
    nC2 <- n*(n-1) / 2
    combn2 <- combn(n,2)
    y <- yMat[t(combn2)]
    if(!isTRUE(sparse)){
        X <- matrix(0,nC2,n)
        X[cbind(1:length(combn2[1,]), combn2[1,])] <- 1L
        X[cbind(1:length(combn2[2,]), combn2[2,])] <- 1L
    }
    list(X = X, y = y)
}


designUnrooted <- function(tree,order=NULL){
    if(is.rooted(tree))tree = unroot(tree)
    p=bipartition(tree)
    if(!is.null(order)) p=p[,order]
    n = dim(p)[1]
    m = dim(p)[2]
    res = matrix(0,(m-1)*m/2, n)
    k = 1
    for(i in 1:(m-1)){
        for(j in (i+1):m){
            res[k,] = p[,i]!=p[,j]
            k=k+1
        }
    }
    colnames(res) = paste(tree$edge[,1],tree$edge[,2],sep="<->")
    res
}


designUltra <- function (tree, sparse=FALSE)
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
        tree = reorderPruning(tree)
    leri = allChildren(tree)
    bp = bip(tree)
    n = length(tree$tip)
    l = tree$Nnode
    nodes = integer(l)
    k = 1L
    u=numeric( n * (n - 1)/2)
    v=numeric( n * (n - 1)/2)
    m = 1L
    for (i in 1:length(leri)) {
        if (!is.null(leri[[i]])) {
            ind =  getIndex(bp[[leri[[i]][1] ]], bp[[leri[[i]][2] ]], n)
            li = length(ind)
            v[m: (m+li-1)]=k
            u[m: (m+li-1)]=ind
            nodes[k]=i
            m = m+li
            k = k + 1L
        }
    }
    if(sparse) X = sparseMatrix(i=u,j=v, x=2L)
    else{
        X = matrix(0L, n * (n - 1)/2, l)
        X[cbind(u,v)]=2L
    }
    colnames(X) = nodes
    attr(X, "nodes") = nodes
    X
}


reorderPruning <- function (x, ...)
{
    parents <- as.integer(x$edge[, 1])
    child <- as.integer(x$edge[, 2])
    root <- as.integer(parents[!match(parents, child, 0)][1])  # unique out
    if (length(root) > 2)
        stop("more than 1 root found")
    n = length(parents)
    m = max(x$edge)  # edge  parents
    neworder = .C("reorder", parents, child, as.integer(n), as.integer(m), integer(n), as.integer(root-1L), DUP=FALSE)[[5]]
    x$edge = x$edge[neworder,]
    x$edge.length = x$edge.length[neworder]
    attr(x, "order") <- "pruningwise"
    x
}
