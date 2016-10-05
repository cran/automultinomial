#' Cholesky decomposition for positive semidefinite matrix
#'
#' Takes a square symmetric positive definite or positive semidefinite matrix and returns a Cholesky factor.
#' @param a A matrix
#' @keywords cholesky, semidefinite
#' @export
#' @return The lower triangular Cholesky factor
#' @examples 
#' X=replicate(100,rnorm(99))
#' 
#' #a is positive semidefinite but not positive definite
#' a_0=t(X)%*%X
#' b=chol_psd(a_0)
#' 
#' #recover a
#' a_1=b%*%t(b)
#' diff=max(abs(a_1-a_0))
chol_psd<-function(a){
  n = dim(a)[1];
  root = matrix(0,n,n);
  
  for (i in 1:n){
    sum = 0;
    if (i>1){
      sum = sum(root[i,1:(i-1)]^2);
    }
    
    x = a[i,i] - sum;
    
    if (x<0){
      x = 0;
    }
    
    root[i,i] = sqrt(x);
    
    if (i < n){
      for (j in (i+1):n){
        
        if (root[i,i] == 0){
          x=0;
        }
        else{
          sum = 0;
          if (i>1) {
            sum = root[i,1:(i-1)] %*% t(t(root[j,1:(i-1)]))
          }
          x = (a[i,j] - sum)/root[i,i];
        }
        
        root[j,i] = x;
      }
    }
  }
  return(root);
}