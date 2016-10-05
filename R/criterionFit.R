#' Minimize criterion (AIC or BIC)
#'
#' Takes a list of coefficient values and finds the values which minimize AIC or BIC
#' @param thetaList a list of coefficient values (probably from a LASSO path)
#' @param X a design matrix
#' @param z response matrix (for multinomial or logistic data)
#' @param A a list of adjacency matrices
#' @param BIC logical, use BIC if true or AIC if false
#' @param logLike function to compute the loglikelihood: either centered or uncentered
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices internal indices
#' @keywords AIC, BIC
#' @return List containing a minimum criterion coefficient vector and the criterion for that vector
#' @examples
#' #criterionFit(thetaList,X,z,A,BIC,logLike,constraint,thetaInds,innerIndices)
criterionFit<-function(thetaList,X,z,A,BIC,logLike,constraint,thetaInds,innerIndices){
  logLik=NULL
  k=dim(z)[2]
  nTheta=dim(thetaList)[1]
  nEta=(k-1)
  nBeta=length(thetaList[[1]])[1]/k-nEta
  bestBIC=10^30
  constant=log(dim(z)[1])
  if (BIC==FALSE){constant=2}
  BICtemp=0
  bestTheta=thetaList[[1]]
  nTheta=length(thetaList)
  BICvec=1:nTheta
  for (i in 1:nTheta){
    dev=-2*logLike(c(thetaList[[i]]),X,z,A,constraint,thetaInds,innerIndices)
    #we don't really have k*p parameters
    temp=matrix(thetaList[[i]],ncol=k)
    nBetaPar=sum(apply(abs(temp[1:nBeta,]),1,sum)!=0)*(k-1)
    if (constraint=="symmetric"){
      nEtaPar=k*(k-1)/2
    }
    if (constraint=="diagonal"){
      nEtaPar=k-1
    }
    nPar=nBetaPar+nEtaPar
    penalty=constant*nPar
    BICtemp=dev+penalty
    BICvec[i]=BICtemp
    if (BICtemp<bestBIC){
      bestTheta=thetaList[[i]]
      bestBIC=BICtemp
    }
  }
  return(list(bestTheta,criterion=bestBIC))
}