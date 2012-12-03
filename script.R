library(phangorn)
library(inline)
library(RcppArmadillo)

## postRcpp <- cxxfunction(signature(data ="list",  P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/likelihood.cpp"), collapse = "\n" ))

## cvRcpp <- cxxfunction(signature(data ="list", cvData="list",P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/cvlikelihood.cpp"), collapse = "\n" ))

cvllobs <- cvCluster(cvDataListobs, 267, maxCluster = 6, mCV=200, beta=0.5, initPList = initPList)
apply( X =cvllobs, 2, FUN = mean, na.rm = TRUE)

resultHR <- list()
for(i in 1:6){
    initDataP <- initPList[[i]]
    resultHR[[i]] <- mtDataResult <- mtCluster(finalPhy, obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat,  initP=initPList[[i]], maxIter = 100, tol = 1e-6)
}

DNAresult <- list()
for(i in 1:6){
    initDataP <- initPList[[i]]
    DNAresult <-  DNACluster(finalPhy,  initP=initPList[[i]], maxIter = 100, tol = 1e-6)
}

save(result, cvllobs, file= "resultHR.RData")

