#' Compute log likelihood
#'
#' Compute log likelihood for uncentered or independent case models
#' @param theta a matrix of coefficinets
#' @param X design matrix
#' @param z response matrix
#' @param A list of adjacency matrices/list containing adjacency matrix
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices internal grid indices
#' @keywords adjacency, spatial, neighborhood, grid
#' @export
#' @return A sparse adjacency matrix
#' @examples 
#' #logLikelihood(theta,X,z,A,constraint,thetaInds)
logLikelihood<-function(theta,X,z,A,constraint,thetaInds,innerIndices){
  k=dim(z[2])
  n=dim(z[1])
  
  k=dim(z)[2]
  n=dim(z[1])
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  
  #enforce constraint on eta's
  if (constraint !="none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  theta=matrix(theta,ncol=k)
  
  #new X matrix
  newX=as.matrix(cbind(X,A[[1]]%*%z[,2:k]))
  
  first=newX%*%theta
  first=sum(first*z[innerIndices,])
  second=exp(newX%*%theta)
  second=apply(second,1,sum)
  second=sum(log(second))
  loglik=first-second
  return(loglik)
}