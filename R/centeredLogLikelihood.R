#' Centered model log pseudolikelihood
#' 
#' Computes a log pseudolikelihood for the centered model, taking coefficients, a design matrix, a response matrix, and a list of adjacency matrices 
#' @param theta a matrix of coefficients with dimension (p+(k-1)*length(A)) by k; frequently A will have length 1
#' @param X nxp design matrix
#' @param z nxk response matrix
#' @param A a list of adjacency matrices
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices internal grid indices
#' @keywords centered, gradient
#' @return a log (pseudo-)likelihood
#' @export
#' @examples 
#' #g=centeredLogLikelihood(theta,X,z,A)
centeredLogLikelihood<-function(theta,X,z,A,constraint,thetaInds,innerIndices){
  k=dim(z)[2]
  n=dim(z[1])
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  
  #enforce constraint on eta's
  if (constraint != "none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  theta=matrix(theta,ncol=k)
  newX=totalX(theta,X,z,A)
  first=newX%*%theta
  first=sum(first*z)
  second=exp(newX%*%theta)
  second=apply(second,1,sum)
  second=sum(log(second))
  loglik=first-second
  return(loglik)
}

