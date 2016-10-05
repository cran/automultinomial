#' Compute sandwich variance estimate
#'
#' Computes the sandwich variance estimate for a given neighborhood structure and coefficient estimate
#' @param theta a coefficient matrix
#' @param X a design matrix
#' @param z a matrix of multinomial responses
#' @param A a list of adjacency matrices
#' @param centered use centered or uncentered model
#' @param I Fisher information matrix
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices internal grid indices
#' @keywords adjacency, spatial, neighborhood, grid
#' @return a covariance matrix
#' @examples 
#' #
sandwichVariance<-function(theta,X,z,A,centered,I=NULL,constraint,thetaInds,innerIndices){
  if (centered){
    J=JMatrixCentered(theta,X,z,A,constraint,thetaInds)
  }
  else{
    J=JMatrixUncentered(theta,X,z,A,constraint,thetaInds,innerIndices)
  }
  var=MASS::ginv(I)%*%J%*%MASS::ginv(I)
  return(var)
}