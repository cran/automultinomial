#' Block coordinate descent
#'
#' Computes, with group lasso regularization penalty, a lasso path for a model with user specified loglikelihood gradient and hessian.
#' @param X design matrix
#' @param z response matrix
#' @param groups a grouping of coefficients
#' @param penaltyFactor a vector specifying the relative penalty level for each coefficient group
#' @param nLambda number of penalty parameters on the lasso path
#' @param H optional Hessian at loglikelihood maximum, used for the quadratic approximation in the centered model
#' @param theta optional coefficients at loglikelihood maximum, used for the quadratic approximation in the centered model
#' @keywords block, coordinate, descent
#' @return list containing the lambdas on the coordinate descent path and the coefficients for each lambda
#' @examples 
#' #
multinomialBCD<-function(X,z,groups,penaltyFactor,nLambda,H=NULL,theta=NULL){

  #convergence tolerance
  tol=10^-7
  #compute starting point for lambda
  #get number of groups
  categories=unique(groups)
  g=length(categories)
  k=dim(z)[2]
  p=dim(X)[2]
  if (!is.null(H)){
    p=dim(H)[2]/k
  }
  #get number of coefficients in each group
  groupLengths=rep(0,g)
  for (i in 1:g){
    groupLengths[i]=sum(groups==categories[i])
  }
  
  #starting beta
  beta=matrix(0,p,k)
  
  #starting derivative norms
  norms=1:sum(penaltyFactor>0)
  count=0
  for (i in 1:g){
    if (penaltyFactor[i]>0){
      count=count+1
      cat=categories[i]
      inds=(groups==cat)
      gradient=-H[inds,]%*%(c(beta-theta))
      norms[count]=sum(gradient^2)^0.5
    }
  }
  lambdas=norms/(groupLengths[penaltyFactor>0]*k*penaltyFactor[penaltyFactor>0])^0.5
  lambdaStart=max(lambdas)
  num=nLambda
  lambdaVec=lambdaStart*(num:1/num)
  betaList=list()
  
  
  #want the negative hessian to minimize -loglik+penalty
  H=-H
  
  #compute list of eigenvalues of the group Hessian matrices
  eigList=1:g
  for (i in 1:g){
    cat=categories[i]
    inds=(groups==cat)
    eigList[i]=max(eigen(H[inds,inds])$values)
  }
  
  oldBeta=beta
  #for each lambda
  count=0
  pb=utils::txtProgressBar(min=1,max=length(lambdaVec),initial=1,width=60,style=3)
  for (lambda in lambdaVec){
    count=count+1
    diff=tol+1
    while(diff>tol){
      diff=0
      #cycle through groups
      for (i in 1:g){
        
        #get coefficient indices for this group
        cat=categories[i]
        inds=(groups==cat)
        hess=H[inds,inds]
        betaK=c(beta[inds])
        gradient=(H[inds,]%*%(c(beta-theta)))-hess%*%betaK
        maxEig=eigList[i]
        betaR=-1/maxEig*(gradient+hess%*%betaK)+betaK
        lambda1=lambda*penaltyFactor[i]
        betaK=max(0,(1-lambda1/(maxEig*sum(betaR^2)^0.5)))*betaR
        beta[inds,]=betaK
      }
      diff=max(abs(beta-oldBeta))
      oldBeta=beta
    }
    #print(count)
    betaList[[count]]=beta
    utils::setTxtProgressBar(pb,count)
  }
  cat("\n")
  return(list(betaList=betaList,lambda=lambdaVec))
}