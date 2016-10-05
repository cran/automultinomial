#' Uncentered gradient
#' 
#' Takes coefficient values and a design matrix and computes the gradient of the loglikelihood
#' @param theta matrix of coefficients
#' @param X an nxp design matrix
#' @param z an nxk response matrix
#' @param A list of adjacency matrices/list containing adjacency matrix
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices the (1-dim.) location/index for internal points in the grid
#' @keywords multinomial, gradient
#' @export
#' @return a gradient vector
#' @examples 
#' #gradient=uncenteredLogLikGrad(theta,X,z,A,constraint,thetaInds)
uncenteredLogLikGrad<-function(theta,X,z,A,constraint,thetaInds,innerIndices){
  X=as.matrix(X)
  if (constraint != "none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  n=dim(z)[1]
  k=dim(z)[2]
  p=length(c(theta))/(k)
  X=cbind(X,A[[1]]%*%z[,2:k])
  X=as.matrix(X)
  grad=rep(0,p*k)
  theta=matrix(theta,ncol=k)
  
  mu=multProb(X,theta)
  grad=rep(0,p*k)
  for (i in 1:k){
    start=(i-1)*p+1
    temp=t(X)%*%(z[innerIndices,i]-mu[,i])
    grad[start:(start+p-1)]=temp
  }
  k=dim(z)[2]
  val=matrix(grad,ncol=k)
  nBeta=dim(X)[2]-(k-1)
  nEta=length(A)*(k-1)
  etaMat=val[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  etaMat[,1]=0
  if (constraint=="none"){
    return(val)
  }
  if (constraint=="diagonal"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=diag(diag(etaSubMat),nrow=dim(etaSubMat)[2],ncol=dim(etaSubMat)[2])
    }
    val[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  if (constraint=="symmetric"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaSubMat[upper.tri(etaSubMat)]=etaSubMat[upper.tri(etaSubMat)]+t(etaSubMat)[upper.tri(t(etaSubMat))]
      etaSubMat[lower.tri(etaSubMat)]=0
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=etaSubMat
    }
    val[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  return(val)
}