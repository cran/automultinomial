#' Logistic autoregression
#' 
#' Fits an autologistic or automultinomial logit model by pseudolikelihood.
#' Fits the traditional uncentered autologistic model by default, or its multicategory analog. 
#' Also capable of fitting the centered autologistic/automultinomial models. Performs an optional group LASSO variable 
#' selection step, with an option to group related variables. For the automultinomial model, the variable selection step
#' treats coefficients relating to a single predictor as grouped, as in glmnet with the option type.multinomial="grouped". This means that
#' entire rows of the multinomial coefficient matrix are selected or removed.
#' 
#' @param X nxp design matrix, where the first column contains an intercept
#' @param z nxk response matrix or a length n vector containing the response type of each observation
#' @param A an adjacency matrix specifying the neighborhood structure
#' @param innerIndices vector containing points to be treated as internal to the grid (defaults to all points if unspecified)
#' @param groups the grouping of the (non-intercept) coefficients. A vector of length (p-1)
#' @param nLambda the length of the lasso path, defaults to 101
#' @param centered logical, use centered model (TRUE) or uncentered model (FALSE)? Defaults to FALSE
#' @param BIC logical, default TRUE uses BIC and FALSE uses AIC for variables selection
#' @param select logical, do variable selection? defaults to TRUE to include a variable selection step
#' @param standardize logical, standardize non-intercept columns? Default is TRUE and usually makes more sense
#' @param constraint the constraint on the eta matrix to ensure a valid joint distribution: "diagonal" (default) or "symmetric" 
#' @param NR use Newton-Raphson optimizer? Defaults to false, but NR may give faster and sharper convergence for uncentered model
#' @keywords logistic, multinomial, regression, autoregression, centered, uncentered, variable, selection, autologistic automultinomial
#' @export
#' @return a list containing model fits and variance estimates 
#' @examples
#' \dontrun{
#' #Simulate data and estimate coefficients for binary data (example 1) 
#' #and multicategory data (example 2)
#' 
#' 
#' ###### Example 1: binary data
#' 
#' m=100
#' 
#' #Generate adjacency matrix for first spatial nearest neighbors on 100x100 grid
#' #with boundary conditions
#' 
#' #the first element of A is an adjacency matrix; the second contains the locations
#' #of internal grid points
#' A=adjMat(m,m,boundary=TRUE)
#' 
#' #generate data and design matrix using example coefficient values
#' set.seed(42)
#' 
#' #some of the predictor values are zero: these variables are irrelevant to the response.
#' #In the model fitting step, setting select=TRUE will 
#' #use the LASSO to attempt to select the relevant variables.
#' beta=cbind(c(-0.1,0.5,-0.5,0.4,-0.4,rep(0,7)))
#' eta=0.5
#' 
#' theta=c(beta,eta)
#' 
#' #X matrix with 12 predictors and 10000 observations for a 100x100 grid
#' X=cbind(rep(1,m^2),replicate(11,rnorm(m^2)))
#' 
#' #generate data for 100x100 grid with "0" boundary conditions
#' z1=simulateData(theta=theta,X=X,m,m,centered=FALSE,k=2,burn=30,draws=1,boundary=0)
#' 
#' #view responses with plotGrid
#' plotGrid(z1,m+2,m+2)
#' 
#' #model fitting
#' model1=MPLE(X,z1,A=A$A,innerIndices=A$innerIndices,nLambda=101,centered=FALSE,BIC=TRUE,select=TRUE)
#' 
#' #significance testing
#' pValues1=significanceTest(model1)
#' 
#' 
#' ###### Example 2: Data with more than two response categories 
#' 
#' #i) simulate data from the automultinomial (Potts) model with 4
#' #response categories on a 100x100 grid without boundary conditions
#' 
#' #ii) model fitting
#' 
#' m=100
#' #Generate adjacency matrix for first spatial nearest neighbors on 100x100 grid
#' A=adjMat(m,m,boundary=FALSE)$A
#' 
#' #generate data and design matrix using example coefficient values
#' set.seed(42)
#' #for multinomial data with k categories, beta and eta will have (k-1) columns.
#' #in this example, some of the rows of beta, corresponding 
#' #to a single variable in the X matrix, are zero. 
#' #These variables are irrelevant to the response. In the model fitting step, 
#' #setting select=TRUE will use the group LASSO to attempt to select the relevant variables.
#' 
#' beta=cbind(c(0.1,0.2,0.4,-0.3,0.5,rep(0,7)),c(-0.5,0.5,-0.5,0.2,-0.4,rep(0,7)),
#' c(-0.3,0.2,0.7,-0.3,0.2,rep(0,7)))

#' #for multicategory data, eta[i,j] corresponds to the effect of neighbors of type i 
#' #on the response probability of type j, treating the reference category as response type 0. 
#' #eta is assumed to be a diagonal or symmetric (k-1) by (k-1) matrix. 
#' eta=cbind(c(0.5,-0.1,-0.3),c(-0.1,0.4,0.1),c(-0.3,0.1,0.5))
#' 
#' theta=rbind(beta,eta)
#' 
#' #X matrix with 12 predictors and 10000 observations for a 100x100 grid
#' X=cbind(rep(1,m^2),replicate(11,rnorm(m^2)))
#' 
#' #generate data for 100x100 grid.
#' z2=simulateData(theta=theta,X=X,m,m,centered=FALSE,k=4,burn=30,draws=1)
#'
#' #view responses with plotGrid
#' plotGrid(z2,m,m)
#' 
#' #model fitting
#' 
#' model2=MPLE(X,z2,A,nLambda=101,centered=FALSE,
#' BIC=TRUE,select=TRUE,standardize=TRUE,constraint="symmetric")
#' 
#' #significance testing
#' pValues2=significanceTest(model2)
#' 
#' }
MPLE<-function(X,z,A=NULL,innerIndices=NULL,groups=2:dim(X)[2],nLambda=101,centered=FALSE,BIC=TRUE,select=TRUE,standardize=TRUE,constraint="diagonal",NR=FALSE){
  X=as.matrix(X)
  if (dim(X)[2]==1){
    select=FALSE
    standardize=FALSE
    X[,1]=rep(1,dim(X)[2])
  }
  k=dim(as.matrix(z))[2]
  if (dim(as.matrix(z))[2]==1){
    k=length(unique(z))
    if (length(unique(z))<2){stop("Must have at least two categories to use MPLE")}
  }
  if (is.null(A)){
    stop("Need to input adjacency matrix")
  }
  if (!is.null(A)){
    if (!is.list(A)){
      A=list(A)  
    }
  }
  if (is.null(innerIndices)){
    innerIndices=1:dim(A[[1]])[1]
  }
  A[[1]]=A[[1]][innerIndices,]


  #adjust variances and coefficients to new parameterization
  #change of parameterization is linear so new variance is just A%*%theta%*%A'
  #we only want to get rid of the variance components of the first p predictors of this new matrix
  
  #transformation matrix
  model=MPLEinternal(X,z,A,groups,nLambda,centered,BIC,select,standardize,constraint,NR,innerIndices)
  p=dim(model$MPLE)[1]
  Ashift=matrix(0,k*p,k*p)
  for (i in 1:(k-1)){
    for (j in 1:p){
      Ashift[i*p+j,j]=-1
      Ashift[i*p+j,i*p+j]=1
    }
  }
  
  #shift MPLE
  MPLEtemp=(Ashift%*%c(model$MPLE))[(p+1):(k*p)]
  MPLEtemp=matrix(MPLEtemp,ncol=k-1)
  
  #shift LASSO path
  thetaList=NULL
  theta=NULL
  theta_variance=NULL
  subModelBIC=NULL
  if (select==TRUE){
    lassoPath=model$lassoPath
    thetaList=list()
    
    #convert lasso path to ordinary reference category parameterization
    for (i in 1:length(lassoPath)){
      thetaList[[i]]=(Ashift%*%c(lassoPath[[i]]))[(p+1):(k*p)]
    }
  }
  #shift variance
  MPLE_BIC=model$MPLE_Criterion
  MPLE_variance=(Ashift%*%model$MPLE_variance%*%t(Ashift))[(p+1):(k*p),(p+1):(k*p)]
  
  #variance of lower triangle of eta equals variance of upper triangle (covariance less than full rank)
  if (constraint=="symmetric"){
    MPLE_variance=varianceTransformation(MPLEtemp,MPLE_variance)
  }
  
  #now deal with the subModels
  if (select==TRUE){
    theta=model$subModel
    theta_variance=model$subModelVariance
    if (!is.null(model$subModel)){
      theta=(Ashift%*%c(theta))[(p+1):(k*p)]
      p=sum(theta!=0)/(k-1)
      if (constraint!="symmetric"){
        p=((sum(theta!=0))-(k-1)+(k-1)*(k-1))/(k-1)
      }
      theta=matrix(theta,ncol=k-1)
      Ashift=matrix(0,k*p,k*p)
      for (i in 1:(k-1)){
        for (j in 1:p){
          Ashift[i*p+j,j]=-1
          Ashift[i*p+j,i*p+j]=1
        }
      }
      theta_variance=(Ashift%*%theta_variance%*%t(Ashift))[(p+1):(k*p),(p+1):(k*p)]
      if (constraint=="symmetric"){
        nonZeroInds=which(apply(abs(theta),1,sum)!=0)
        theta_variance=varianceTransformation(theta[nonZeroInds,],theta_variance)
      }
    }
    subModelBIC=model$subModelCriterion
  }
  return(list(lassoPath=thetaList,MPLE=MPLEtemp,MPLE_variance=as.matrix(MPLE_variance),MPLE_Criterion=MPLE_BIC,subModel=theta,subModelVariance=theta_variance,subModelCriterion=subModelBIC))
}