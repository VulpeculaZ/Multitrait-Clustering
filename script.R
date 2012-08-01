library(phangorn)
library(inline)
library(RcppArmadillo)

## postRcpp <- cxxfunction(signature(data ="list",  P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/likelihood.cpp"), collapse = "\n" ))

## cvRcpp <- cxxfunction(signature(data ="list", cvData="list",P="matrix", contrast="matrix", nrs="integer" , ncs="integer", ncos="integer", bfs="numeric", ecps="matrix"), plugin="RcppArmadillo", body=paste( readLines("./Source/cvlikelihood.cpp"), collapse = "\n" ))

## cvllobs <- cvCluster(cvDataListobs, 267, maxCluster = c(1:7), mCV=100, beta=0.5, initPList = initPList)

resultHR <- list()
for(i in 1:7){
    initDataP <- randomPMat(0.05, i, 267)
    resultHR[[i]] <- hostCluster(obsHost = hostDistMat$obsHost, distMat=hostDistMat$distMat, initP=initDataP, maxIter = 200, tol = 1e-4)
}

save(result, cvllobs, file= "resultHR.RData")

