#' Multinomial regression
#' 
#' Fits an automultinomial model (by pseudolikelihood). Fits the uncentered auto-model by default, 
#' #when an adjacency matrix is given. Performs an optional variable selection step, treating each
#' #variable as a group for all response categories, with an option to group related variables.
#' @param X nxp design matrix, where the first column must contain an intercept
#' @param z nxk response indicator matrix, or a length n vector with the response type of each observation
#' @param A an adjacency matrix or list of adjacency matrices defining neighborhood structures
#' @param groups the grouping of the (non-intercept) coefficients. A vector of length (p-1)
#' @param nLambda the length of the lasso path, defaults to 101
#' @param centered logical, TRUE for centered and FALSE for uncentered, default is TRUE
#' @param BIC logical, default TRUE uses BIC and FALSE uses AIC
#' @param select logical, TRUE for variable selection step, default is TRUE
#' @param standardize logical, standarize non-intercept columns? Default is TRUE, FALSE is almost certainly a bad idea
#' @param constraint "symmetric" or "diagonal"
#' @param polish use full Newton-Raphson to sharpen convergence from BFGS
#' @param innerIndices indices of internal grid points
#' @keywords multinomial, regression, autoregression, centered, uncentered, variable, selection
#' @return a list containing the (group) lasso path, 
#' @examples
#' #
MPLEinternal<-function(X,z,A=NULL,groups=2:dim(X)[2],nLambda=101,centered=TRUE,BIC=TRUE,select=TRUE,standardize=TRUE,constraint,polish,innerIndices){

  #standardize the X matrix
  X=as.matrix(X)
  
  #intercept only model
  if (dim(X)[2]==1){
    X[,1]=rep(1,dim(X)[1])
    standardize=FALSE
    groups=NULL
    select=FALSE
  }
  
  #save unstandardized X matrix
  tempX=X
  
  #use AIC or BIC?
  string="BIC"
  if (BIC==FALSE){string="AIC"}
  
  #turn response vector into matrix
  if (dim(as.matrix(z))[2]==1){
    z=factor(z)
    z=as.numeric(z)
    k=max(z)
    n=length(z)
    temp=matrix(0,n,k)
    for (i in 1:n){
      temp[i,z[i]]=1
    }
    z=temp
  }
  
  #get number of categories and observations
  k=dim(z)[2]
  n=dim(z)[1]
  
  #standardizing, making sure first column is intercept
  if (dim(X)[2]>=2){
    means=apply(X[,2:dim(X)[2],drop=FALSE],2,mean)
    sds=apply(X[,2:dim(X)[2],drop=FALSE],2,stats::sd)
    X[,1]=rep(1,dim(X)[1])
    if (standardize){
      X[,2:dim(X)[2]]=scale(X[,2:dim(X)[2]])
    }
  }

  #get nBeta and nEta
  nBeta=dim(X)[2]
  #only (k-1)*length(A) identifiable etas
  nEta=length(A)*(k-1)
  theta=matrix(0,nBeta+nEta,k)
  
  #get indices for later use
  thetaInds=indsFunction(theta,X,z,A)
  
  cat("Fitting Maximum Pseudolikelihood Estimate\n")
  #fitting with standardized X matrix
  logLike=NULL
  gr=NULL
  hess=NULL
  if (centered){
    logLike=centeredLogLikelihood
    gr=centeredLogLikGrad
    hess=centeredLogLikHess
  }
  if (!centered){
    logLike=logLikelihood
    gr=uncenteredLogLikGrad
    hess=uncenteredLogLikHess
  }
  
  #BFGS first to get into pos. def. quadratic region for centered case
  result=NULL
  if ((!polish | centered)){
    result=stats::optim(c(theta),logLike,gr=gr,X,z,A,constraint,thetaInds,innerIndices,method=c("BFGS"),control=list(fnscale=-1,maxit=2000))  
    theta=toConstraint(result$par,X,z,A,constraint)
    theta=matrix(theta,ncol=k)
  }
  
  ###################
  iterations=1
  if (polish){
    maxGrad=1
    while ((maxGrad/(dim(X)[1]))>10^-11 & iterations<=100){
      #print(logLike(theta,X,z,A,constraint,thetaInds))
      iterations=iterations+1
      H=hess(theta,X,z,A,constraint,thetaInds)
      gradient=gr(c(theta),X,z,A,constraint,thetaInds,innerIndices)
      #print(matrix(gradient,ncol=k))
      step=(MASS::ginv(H))%*%c(gradient)
      theta=c(theta)
      theta=theta-step
      maxGrad=max(abs(gradient))
      #print(maxGrad)
      theta=toConstraint(theta,X,z,A,constraint)
    }
    theta=matrix(theta,ncol=k)
  }
  ###################
  activeInds=NULL
  if (constraint=="symmetric"){
    activeInds=thetaInds$betaInds+thetaInds$etaDiag+thetaInds$etaOffDiagUpper
  }
  if (constraint=="diagonal"){
    activeInds=thetaInds$betaInds+thetaInds$etaDiag
  }
  activeInds=as.logical(activeInds)
  cat("Checking for convergence\n")
  H=hess(theta,X,z,A,constraint,thetaInds)
  maxEigen=eigen(H[activeInds,activeInds])$values[nBeta+1]
  
  convergence=FALSE
  if (!polish){
    convergence=(result$convergence==0)
  }
  if ((!polish & !convergence) | (polish & iterations>100)){stop("Optimization failed, sample size is probably too small for centered model.")}
  if (maxEigen>=0){stop("Optimization failed, sample size is probably too small for centered model.")}
  MPLE=matrix(theta,ncol=k)
  #MPLE=MPLE-apply(MPLE,1,mean)
  
  #compute criterion for MPLE
  MPLEdev=-2*logLike(theta,X,z,A,constraint,thetaInds,innerIndices)
  constant=log(n)
  #use AIC if not BIC
  if (!BIC){constant=2}
  #only (k-1) columns of identifiable parameters
  MPLE_BIC=MPLEdev+constant*(nBeta*(k-1)+k*(k-1)/2*(constraint=="symmetric")+(k-1)*(constraint=="diagonal"))
  
  #compute variance
  cat("Computing MPLE variance\n")
  
  #compute variance on original scale, MPLEtemp will be on original scale and MPLE is on standardized scale
  MPLEtemp=MPLE
  #convert coefficients to unstandardized X matrix coefficients
  if (standardize){
    MPLEtemp=unStandardize(list(MPLE),means,sds,nRow=dim(X)[2])[[1]]
  }
  
  #use unstandardized X matrix tempX
  #MPLE_variance=NULL
  #Fisher information on original scale
  Htemp=hess(MPLEtemp,tempX,z,A,constraint,thetaInds)
  MPLE_variance=sandwichVariance(MPLEtemp,tempX,z,A,centered=centered,I=-Htemp,constraint,thetaInds,innerIndices)
  
  ####################
  theta=NULL
  theta_variance=NULL
  fullTheta=list(NULL,NULL)
  thetaList=NULL
  subModelBIC=NULL
  
  #if variable selection should be done
  if (select==TRUE){  
    #LASSO step
    #BCD
    #get magnitudes of each group
    groups=as.numeric(factor(groups))
    groups=c(0,groups,rep(0,nEta))
    
    cats=unique(groups)
    penaltyLevels=rep(0,length(cats))
    
    #skip intercept
    count=0
    #adaptive lasso
    for (i in cats){
      count=count+1
      inds=(groups==i)
      coefs=c(MPLE[inds,])
      magnitude=sum(coefs^2)^0.5
      #don't penalize intercept
      penaltyLevels[count]=(count!=1)*1/magnitude
    }
    
    #LASSO step
    cat("LASSO variable selection step\n")
    penaltyFactor=penaltyLevels
    lassoObj=multinomialBCD(X,z,groups,penaltyFactor,nLambda=nLambda-1,H,MPLE-MPLE[,1])
    
    #the lasso only optimizes one part of the eta matrix in the symmetric parameterization
    #the Hessian and gradient used in the lasso step takes the other eta component into account
    lambda=lassoObj$lambda
    thetaList=lassoObj$betaList
    thetaList[[nLambda]]=MPLE-MPLE[,1]
    for (i in 1:length(thetaList)){
      thetaList[[i]]=toConstraint(thetaList[[i]],X,z,A,constraint)
      thetaList[[i]]=matrix(thetaList[[i]],ncol=k)
    }
    
    #criterion minimization
    fullTheta=criterionFit(thetaList,X,z,A,BIC,logLike,constraint,thetaInds,innerIndices)
    theta=matrix(fullTheta[[1]],ncol=k)
    #get nonzero indices
    thetaNonZero=apply(abs(theta),1,sum)!=0
    inds=apply(abs(theta[1:dim(X)[2],]),1,sum)!=0
    nBeta=sum(inds)
    #refit on original scale
    cat("Refitting\n")
    
    #scale lasso'd parameters back
    if (standardize){
      theta=unStandardize(list(theta),means,sds,nRow=dim(X)[2])[[1]]
      thetaList=unStandardize(thetaList,means,sds,nRow=dim(X)[2])
    }
    thetaInds=indsFunction(theta[thetaNonZero,,drop=FALSE],tempX[,inds,drop=FALSE],z,A)
    if (!polish){
      result=stats::optim(c(theta[thetaNonZero,,drop=FALSE]),logLike,gr=gr,tempX[,inds,drop=FALSE],z,A,constraint,thetaInds,innerIndices,method=c("BFGS"),control=list(fnscale=-1,maxit=2000,reltol=10^-12))  
      theta[thetaNonZero,]=toConstraint(result$par,X[,inds,drop=FALSE],z,A,constraint)
    }
    theta=matrix(theta,ncol=k)
    
    ###################
    iterations=1
    if (polish){
      maxGrad=1
      while ((maxGrad/(dim(X)[1]))>10^-11 & iterations<=100){
        H=hess(theta[thetaNonZero,,drop=FALSE],tempX[,thetaNonZero[1:(length(thetaNonZero)-(k-1))],drop=FALSE],z,A,constraint,thetaInds)
        gradient=gr(c(theta[thetaNonZero,,drop=FALSE]),tempX[,thetaNonZero[1:(length(thetaNonZero)-(k-1))],drop=FALSE],z,A,constraint,thetaInds,innerIndices)
        step=(MASS::ginv(H))%*%c(gradient)
        start=c(theta[thetaNonZero,])
        theta[thetaNonZero,]=start-step
        maxGrad=max(abs(gradient))
        theta=toConstraint(theta,X,z,A,constraint)
        theta=matrix(theta,ncol=k)
        iterations=iterations+1
      }
    }
    ###################
    
    theta=theta-theta[,1]
    cat("Checking for convergence\n")
    H=hess(theta[thetaNonZero,,drop=FALSE],tempX[,inds,drop=FALSE],z,A,constraint,thetaInds)
    activeInds=NULL
    if (constraint=="symmetric"){
      activeInds=thetaInds$betaInds+thetaInds$etaDiag+thetaInds$etaOffDiagUpper
    }
    if (constraint=="diagonal"){
      activeInds=thetaInds$betaInds+thetaInds$etaDiag
    }
    activeInds=as.logical(activeInds)
    H1=H[activeInds,activeInds]
    maxEigen=eigen(H1)$values[nBeta+1]
    convergence=FALSE
    if (!polish){
      convergence=(result$convergence==0)
    }
    if ((!polish & !convergence) | (polish & iterations>100)){"Refitting failed, using partially optimized fit from lasso path"}
    
    #compute BIC
    nEtaPar=k*(k-1)/2*(constraint=="symmetric")+(k-1)*(constraint=="diagonal")
    nBetaPar=nBeta*(k-1)
    nPar=nEtaPar+nBetaPar
    subModelBIC=-2*logLike(theta,X,z,A,constraint,thetaInds,innerIndices)+constant*nPar
    #compute variance
    cat("Computing submodel variance\n")
    theta_variance=as.matrix(sandwichVariance(theta[thetaNonZero,,drop=FALSE],tempX[,inds,drop=FALSE],z,A,centered=centered,I=-H,constraint,thetaInds,innerIndices=innerIndices))
  }
  return(list(lassoPath=thetaList,MPLE=MPLEtemp,MPLE_variance=as.matrix(MPLE_variance),MPLE_Criterion=MPLE_BIC,subModel=theta,subModelVariance=theta_variance,subModelCriterion=subModelBIC))
}
