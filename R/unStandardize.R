#' Unstandardize coefficients
#'
#' Converts coefficients on a standardized scale to the original scale
#' @param thetaList a list of coefficient matrices
#' @param means means of the non-intercept columns on the original scale
#' @param sds standard deviations of the coefficients on the original scale
#' @param nRow the number of coefficient rows that are affected by standardization and unstandardization
#' @keywords standardize
#' @return a list of coefficient matrices on the original scale
#' @examples 
#' #
unStandardize<-function(thetaList,means,sds,nRow){
  k=dim(thetaList[[1]])[2]
  thetaNew=thetaList
  diff=means/sds
  for (j in 1:length(thetaList)){
    thetaTemp=thetaList[[j]]
    for (i in 1:k){
      thetaNew[[j]][1,i]=thetaTemp[1,i]-sum(thetaTemp[2:nRow,i]*diff)
    }
    thetaNew[[j]][2:nRow,]=thetaTemp[2:nRow,]/sds
  }
  return(thetaNew)
}