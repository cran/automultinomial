#' Add columns to X matrix
#'
#' Adds the A%*%(z-mu) columns to the design matrix. Used for fitting the centered model
#' @param theta matrix of coefficients
#' @param X design matrix
#' @param z matrix of multinomial responses
#' @param A list of adjacency matrices
#' @keywords centered
#' @return an expanded design matrix
#' @examples 
#' #
totalX<-function(theta,X,z,A){
  k=dim(z)[2]
  nEta=length(A)*(k-1)
  nBeta=dim(theta)[1]-nEta
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])[,2:k]
  newX=X
  for (i in 1:length(A)){
    newX=cbind(newX,A[[i]]%*%(z[,2:k]-mu))
  }
  return(as.matrix(newX))
}