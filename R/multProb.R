#' Multinomial probabilities
#' 
#' Takes a design matrix and coefficients and computes the fitted probabilities for these coefficients.
#' @param X nxp design matrix
#' @param beta a pxk coefficient matrix
#' @keywords multinomial, fitted, probabilities
#' @return an nxk matrix of probabilities
#' @examples 
#' #probs=multProb(X,beta)
multProb<-function(X,beta){
  k=dim(beta)[2]
  n=dim(X)[1]
  denominator=rep(0,n)
  for (i in 1:k){
    denominator=denominator+as.vector(exp(X%*%beta[,i]))
  }
  numerator=matrix(0,n,k)
  for (i in 1:k){
    numerator[,i]=as.vector(exp(X%*%beta[,i]))
  }
  prob=numerator/denominator
  return(prob)
}