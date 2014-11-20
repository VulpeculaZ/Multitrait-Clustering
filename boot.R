idenVec <- function(x){
    all(x==x[1])
}
nIdenVec <- function(x,y){
    sum(unlist(x) !=y)
}

bootTheta <- function(theta, y, p, obs, Ez, n, percentile = seq(0, 1, 0.01)){
    yOrig <- y
    y = y - 1L
    logTheta <- log(theta)
    log13Theta <- log(1 - 3 * theta)
    logP <- log(p)
    nY <- dim(y)[1]
    idenY <- apply(X = y, 1, FUN = idenVec)
    diffY <- y[!idenY,]
    lY <- dim(y)[1]
    ldiffY <- dim(diffY)[1]
    M <- length(p)
    simY <- sample.int(M, size = n, replace=TRUE, prob = p)
    boot <- matrix(, n, M)
    nDiff  <- matrix(, nrow=n, ncol=M)
    ## obs Fisher information
    obsDiff <- matrix(, dim(Ez)[1], dim(Ez)[2])
    for(i in 1:M){
        obsDiff[,i] <- mapply(nIdenVec, x=obs, MoreArgs = list(y = yOrig[,i]))
    }
    plogObs <- Ez * (sweep(obsDiff, 2, theta, FUN="/") - sweep(3*(nY-obsDiff), 2, (1-3*theta), FUN="/"))
    csPlogObs <- colSums(plogObs)
    ## Ix|y
    plogObs2 <- csPlogObs %*% t(csPlogObs) - t(plogObs) %*% plogObs
    diag(plogObs2) <- colSums(Ez * (sweep(obsDiff, 2, theta, FUN="/") - sweep(3*(nY-obsDiff), 2, (1-3*theta), FUN="/"))^2)
    ## Ic
    Ic <- diag(colSums(Ez * ( -sweep(obsDiff, 2, theta^2, FUN="/") - sweep(9*(nY-obsDiff), 2, (1-3*theta)^2, FUN="/"))))
    varScore <- plogObs2 - Ic
    ## Bootstrap cutoff
    bootP <- 0
    for(i in 1:M){
        sizeX <- sum(simY == i)
        xSame <- rbinom(n = sizeX, size = lY - ldiffY, prob = 3 * theta[i])
        sampleDiffi <- (sample(c(0,1,2,3), size = sizeX * ldiffY, replace = TRUE, prob = c(1-3*theta[i], rep(theta[i], 3))) + rep(diffY[,i], sizeX)) %% 4
        for(j in 1:M){
            xDiffj <- matrix(sampleDiffi, nrow = ldiffY) != matrix(rep(diffY[,j], sizeX), nrow = ldiffY)
            nDiff[(bootP+1):(bootP + sizeX),j] <- colSums(xDiffj)
        }
        nDiff[(bootP+1):(bootP + sizeX),] <- sweep(nDiff[(bootP+1):(bootP + sizeX),], 1, xSame, FUN = "+")
        bootP <- bootP + sizeX
    }
    logLikeCl <- sweep(nDiff, 2, logTheta, FUN = "*") +  sweep(nY - nDiff, 2, log13Theta, FUN = "*")
    logLikeCl <- sweep(logLikeCl, 2, logP, FUN = "+")
    logLike <- apply(logLikeCl, 1, logsumexp)
    logwj <- sweep(logLikeCl, 1, logLike, FUN = "-")
    boot <- exp(logwj) * sweep(sweep(nDiff, 2, 3* theta * nY, FUN="-"), 2, theta * (1-3 *theta), FUN = "/")
    bootQuad <- rowSums(boot %*% solve(varScore) * boot)
    ## bootMin <- apply(boot, 1, min)
    ## qMin <- quantile(bootMin, probs = percentile)
    qQuad <- quantile(bootQuad, probs = percentile)
    list(qQuad = qQuad, varScore = varScore)
}

bootX <- function(theta, y, p, x, varScore){
    logTheta <- log(theta)
    log13Theta <- log(1 - 3 * theta)
    M <- dim(y)[2]
    nY <- dim(y)[1]
    xDiff <- matrix(, length(x), M)
    for(i in 1:M){
        xDiff[,i] <- mapply(nIdenVec, x=x, MoreArgs = list(y = y[,i]))
    }
    logLikeCl <- sweep(xDiff, 2, logTheta, FUN = "*") + sweep(nY - xDiff, 2, log13Theta, FUN = "*")
    logLikeCl <- sweep(logLikeCl, 2, logP, FUN = "+")
    logLike <- apply(logLikeCl, 1, logsumexp)
    logwj <- sweep(logLikeCl, 1, logLike, FUN = "-")
    boot <- exp(logwj) * sweep(sweep(xDiff, 2, 3* theta * nY, FUN="-"), 2, theta * (1-3 *theta), FUN = "/")
    ## bootMin <- apply(boot, 1, min)
    bootQuad <- rowSums(boot %*% solve(varScore) * boot)
    list(bootQuad = bootQuad)
}

