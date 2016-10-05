#' Uncentered hessian
#' 
#' Takes coefficient values and a design matrix and computes the hessian matrix of the loglikelihood
#' @param theta matrix of coefficients
#' @param X an nxp design matrix
#' @param z response matrix
#' @param A a list of adjacency matrices/list containing an adjacency matrix
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @keywords multinomial, regression
#' @export
#' @return a hessian matrix
#' @examples 
#' #hess=uncenteredLogLikHess(theta,X,z,A)
uncenteredLogLikHess<-function(theta,X,z,A,constraint,thetaInds){
  if (constraint != "none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  X=as.matrix(X)
  k=dim(z)[2]
  p=length(c(theta))/(k)
  n=dim(z)[1]
  fullX=cbind(X,A[[1]]%*%z[,2:k])
  fullX=as.matrix(fullX)
  theta=matrix(theta,ncol=k)
  mu=multProb(fullX,theta)
  hess=matrix(0,p*k,p*k)
  
  for (i in 1:k){
    start1=(i-1)*p+1
    end1=i*p
    
    for (j in 1:k){
      start2=(j-1)*p+1
      end2=j*p
      temp=diag(p)
      if (i==j){
        w1=mu[,i]-mu[,i]^2
        temp=crossprod(fullX,-w1*fullX)
      }
      else{
        w1=-mu[,j]*mu[,i]
        temp=crossprod(fullX,-w1*fullX)
      }
      hess[start1:end1,start2:end2]=temp
    }
  }
  if (constraint=="none"){
    return(hess)
  }
  if (constraint=="symmetric" | constraint=="diagonal"){
    etaHess=hess[thetaInds$etaInds,thetaInds$etaInds,drop=FALSE]
    transform=toConstraintDerivative(theta,X,z,A,constraint)
    etaHess=transform%*%etaHess%*%t(transform)
    betaEtaHess=hess[thetaInds$etaInds,thetaInds$betaInds,drop=FALSE]
    betaEtaHess=transform%*%betaEtaHess
    hess[thetaInds$etaInds,thetaInds$etaInds]=etaHess
    hess[thetaInds$etaInds,thetaInds$betaInds]=betaEtaHess
    hess[thetaInds$betaInds,thetaInds$etaInds]=t(betaEtaHess)
  }
  return(as.matrix(Matrix::forceSymmetric(hess)))
}