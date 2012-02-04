##' <description>
##' Fitting the cluster using non-negative least square.
##' <details>
##' @title cluster.nnls
##' @param y observed value
##' @param X Design matrix
##' @param w Weight
##' @param sparse Is it a sparse matrix?
##' @return Fitted value
##' @author Ziqian Zhou
cluster.nnls <- function(y, X, w=NULL, sparse = FALSE){
    ## need to use nls
    ## nnls() in R base
    ## or nnls for sparse matrices
    ## na.action
    if(any(is.na(y))){
        ind = which(is.na(y))
        X = X[-ind,,drop=FALSE]
        y= y[-ind]
    }
    ## LS solution
    if(sparse == FALSE){
        if(is.null(w))  {
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

combp <- function(X, p){
    return(apply(X, 1, function(x) prod(p[as.logical(x[1:length(p)])])))
}

##' <description>
##' Create the design matrix and format observed values for least square criteria.
##' The dimension of the design matrix will be: n(n-1)(m^2)/2 by nm+m.
##' Observed
##' <details>
##' @title designObsCluster
##' @param n number of observation
##' @param m number of clusters
##' @param yMat Input observed distances. Need to be in the form of
##' distance matrix.
##' @param sparse Will the returned matrix be sparse?
##' @return A design Matrix representing the cluster by a tree structure. And a vector of observed distance.
##' @author Ziqian Zhou
designObsCluster <- function(n, m, yMat, sparse=FALSE){
    nC2 <- n*(n-1) / 2
    combn2 <- combn(n,2)
    y <- yMat[t(combn2)]
    y <- rep(y, m*m)
    if(!sparse){
        X <- matrix(0, n*(n-1)*m^2/2, n*m+m)
        matnC21 <- matnC22<- matnC2 <- matrix(0,nC2,n)
        matnC21[cbind(1:length(combn2[1,]), combn2[1,])] <- 1L
        matnC22[cbind(1:length(combn2[2,]), combn2[2,])] <- 1L
        matnC2 <- matnC21 + matnC22

        for(i in 1:m){
            for(j in 1:m){
                if(i!=j){
                    X[(((i-1)*m + j - 1) * nC2 + 1) : (((i-1)*m + j) * nC2) , ((i-1) * n + 1):(i*n)] <- matnC21
                    X[(((i-1)*m + j - 1) * nC2 + 1) : (((i-1)*m + j) * nC2) , ((j-1) * n + 1):(j*n)] <- matnC22
                    X[(((i-1)*m + j - 1) * nC2 + 1) : (((i-1)*m + j) * nC2) , m*n+i] <- X[(((i-1)*m + j - 1) * nC2 + 1) : (((i-1)*m + j) * nC2) , m*n+j] <- 1
                }
                else{
                    X[(((i-1)*m + j - 1) * nC2 + 1) : (((i-1)*m + j) * nC2) , ((i-1) * n + 1):(i*n)] <- matnC2
                }
            }
        }
        if(m == 2) X <- X[,-(n*m + m)]
    }
    ## else{
    ## }
    list(X = X, y=y)
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
