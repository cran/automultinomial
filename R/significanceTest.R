#' Significance tests for autologistic/automultinomial models
#' 
#' Performs two-sided hypothesis tests for each coefficient in a centered or uncentered
#' autologistic model. The tests use the asymptotic normality of pseudolikelihood estimates 
#' with the Godambe sandwich variance estimator. A test for the spatial correlation eta is also
#' performed. 
#' 
#' For multinomial models, tests are performed for each row of the coefficient matrix, 
#' as each row corresponds to a single predictor. The matrix of eta coefficients is tested 
#' against the null hypothesis of no spatial correlation.
#' @param model a centered or uncentered autologistic fit
#' @export
#' @return a list of p-values and z-values
#' @examples
#' #model=MPLE(...)
#' #significanceTest(model)
significanceTest<-function(model){
  MPLE_pValues=rowPvalues(model$MPLE,model$MPLE_variance)
  subModel_pValues=NULL
  if (!is.null(model$subModel)){
    inds=1:dim(model$subModel)[1]
    inds=inds[apply(abs(model$subModel),1,sum)>0]
    pValues=rowPvalues(model$subModel[inds,,drop=FALSE],model$subModelVariance)
    thetaDims=dim(model$MPLE)
    nBeta=thetaDims[1]-(thetaDims[2])
    nEta=thetaDims[2]
    subModel_betaPvalues=rep(1,nBeta)
    
    betaInds=inds[1:(length(inds)-nEta)]
    subModel_betaPvalues[betaInds]=pValues$beta_pValues
    subModel_pValues=list(beta_pValues=subModel_betaPvalues,eta_pValue=pValues$eta_pValue)
  }
  return(list(MPLE_pValues=MPLE_pValues,subModel_pValues=subModel_pValues))
}


rowPvalues<-function(coefficients,covariance){
  nRow=dim(coefficients)[1]
  nCol=dim(coefficients)[2]
  k=nCol+1
  nEta=k-1
  nBetaRow=nRow-(k-1)
  betaRows=1:nBetaRow
  etaRows=1:nEta+nBetaRow
  
  beta_pValues=1:nBetaRow
  indsMat=matrix(1:(nRow*nCol),nRow,nCol)
  for (i in 1:nBetaRow){
    inds=c(indsMat[i,])
    coefs=c(coefficients[i,])
    covInv=solve(covariance[inds,inds])
    chiSq=t(coefs)%*%covInv%*%coefs
    beta_pValues[i]=1-stats::pchisq(chiSq,nCol)
  }
  
  etaInds=c(indsMat[etaRows,])
  eta_cov=covariance[etaInds,etaInds]
  eta_chiSq=t(coefficients[etaInds])%*%MASS::ginv(eta_cov)%*%coefficients[etaInds]
  eta_pValue=1-stats::pchisq(eta_chiSq,qr(eta_cov)$rank)
  return(pValues=list(beta_pValues=beta_pValues,eta_pValue=eta_pValue))
}