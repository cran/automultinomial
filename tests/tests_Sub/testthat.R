library(testthat)
library(automultinomial)
library(numDeriv)
context("automultinomial")
#Generate adjacency matrices for first and second spatial nearest neighbors on 30x30 grid
n1=adjMat(30,30)

#Convert to list
A=list(n1$A)

#generate data and design matrix using example coefficient values
set.seed(42)
beta=cbind(c(-1,0.75,0.5,1,-0.4,rep(0,7)),c(0.3,-0.25,-3,0,0,rep(0,7)))

#eta
eta=cbind(c(0.5,-0.2),c(-0.2,0.5))

theta=rbind(beta,eta)

#X matrix with 12 predictors and 900 observations for a 30x30 grid
X=cbind(rep(1,900),replicate(11,rnorm(900,0,runif(1,1,2))))

#generate data for 30x30 grid.
zCentered=simulateData(theta,X,30,30,centered=TRUE,k=3,burn=10)
zUncentered=simulateData(theta,X,30,30,centered=FALSE,k=3,burn=10)

theta=theta=cbind(rep(0,14),theta)
thetaInds=indsFunction(theta,X,zCentered,A)

#gradients
test_that("Centered gradient matches numerical gradient", {
  
  t1=grad(centeredLogLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="none",thetaInds)
  t2=c(centeredLogLikGrad(theta,X,zCentered,A,constraint="none",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  t1=grad(centeredLogLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="symmetric",thetaInds)
  t2=c(centeredLogLikGrad(theta,X,zCentered,A,constraint="symmetric",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  t1=grad(centeredLogLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="diagonal",thetaInds)
  t2=c(centeredLogLikGrad(theta,X,zCentered,A,constraint="diagonal",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
})

test_that("Centered hessian matches numerical hessian",{
  t1=jacobian(centeredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="none",thetaInds)
  t2=centeredLogLikHess(theta,X,zCentered,A,constraint="none",thetaInds)
  expect_equal(t1,t2,tolerance=1e-6,check.attributes=FALSE)
  
  t1=jacobian(centeredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="symmetric",thetaInds)
  t2=centeredLogLikHess(theta,X,zCentered,A,constraint="symmetric",thetaInds)
  expect_equal(t1,t2,tolerance=1e-6,check.attributes=FALSE)
  
  t1=jacobian(centeredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zCentered,A,constraint="diagonal",thetaInds)
  t2=centeredLogLikHess(theta,X,zCentered,A,constraint="diagonal",thetaInds)
  expect_equal(t1,t2,tolerance=1e-6,check.attributes=FALSE)
})

test_that("Uncentered gradient matches numerical gradient",{
  t1=grad(logLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,constraint="none",thetaInds)
  t2=c(uncenteredLogLikGrad(theta,X,zUncentered,A,"none",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  t1=grad(logLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,constraint="symmetric",thetaInds)
  t2=c(uncenteredLogLikGrad(theta,X,zUncentered,A,"symmetric",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  t1=grad(logLikelihood,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,constraint="diagonal",thetaInds)
  t2=c(uncenteredLogLikGrad(theta,X,zUncentered,A,"diagonal",thetaInds))
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
})

test_that("Uncentered hessian matches numerical hessian",{
  t1=jacobian(uncenteredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,"none",thetaInds)
  t2=uncenteredLogLikHess(theta,X,zUncentered,A,"none",thetaInds)
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  
  t1=jacobian(uncenteredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,"symmetric",thetaInds)
  t2=uncenteredLogLikHess(theta,X,zUncentered,A,"symmetric",thetaInds)
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
  
  
  t1=jacobian(uncenteredLogLikGrad,theta,method="Richardson",side=NULL,method.args=list(),X,zUncentered,A,"diagonal",thetaInds)
  t2=uncenteredLogLikHess(theta,X,zUncentered,A,"diagonal",thetaInds)
  expect_equal(t1,t2,tolerance=10^-6,check.attributes=FALSE)
})

test_that("Model fits are reasonable for centered case (gradients close to zero)",{
  #fit model with symmetry constraint
  modelCentered=MPLE(X,zCentered,A,nLambda=101,centered=TRUE,BIC=TRUE,select=TRUE,constraint="symmetric",NR=TRUE)
  
  t1=centeredLogLikGrad(cbind(rep(0,14),modelCentered$MPLE),X,zCentered,A,constraint="symmetric",thetaInds)
  subModel=modelCentered$subModel
  
  inds=c(rep(0,14),modelCentered$subModel)!=0
  t2=centeredLogLikGrad(cbind(rep(0,14),modelCentered$subModel),X,zCentered,A,constraint="symmetric",thetaInds)
  
  #hopefully near 0
  expect_equal(c(t1),rep(0,length(t1)),tolerance=10^-6,check.attributes=FALSE)
  expect_equal(c(t2[inds]),rep(0,sum(inds)),tolerance=10^-6,check.attributes=FALSE)
  
  #fit model with diagonal constraint
  modelCentered=MPLE(X,zCentered,A,nLambda=101,centered=TRUE,BIC=TRUE,select=TRUE,constraint="diagonal",NR=TRUE)
  
  t1=centeredLogLikGrad(cbind(rep(0,14),modelCentered$MPLE),X,zCentered,A,constraint="diagonal",thetaInds)
  subModel=modelCentered$subModel
  
  inds=c(rep(0,14),modelCentered$subModel)!=0
  t2=centeredLogLikGrad(cbind(rep(0,14),modelCentered$subModel),X,zCentered,A,constraint="diagonal",thetaInds)
  
  #hopefully near 0
  expect_equal(c(t1),rep(0,length(t1)),tolerance=10^-6,check.attributes=FALSE)
  expect_equal(c(t2[inds]),rep(0,sum(inds)),tolerance=10^-6,check.attributes=FALSE)
})

test_that("Model fits are reasonable for uncentered case (gradients close to zero)",{

  #fit model with symmetry constraint
  modelUncentered=MPLE(X,zUncentered,A,nLambda=101,centered=FALSE,BIC=TRUE,select=TRUE,constraint="symmetric",NR=TRUE)
  
  t1=uncenteredLogLikGrad(cbind(rep(0,14),modelUncentered$MPLE),X,zUncentered,A,constraint="symmetric",thetaInds)
  subModel=modelUncentered$subModel
  
  inds=c(rep(0,14),modelUncentered$subModel)!=0
  t2=uncenteredLogLikGrad(cbind(rep(0,14),modelUncentered$subModel),X,zUncentered,A,constraint="symmetric",thetaInds)
  
  #hopefully near 0
  expect_equal(c(t1),rep(0,length(t1)),tolerance=10^-6,check.attributes=FALSE)
  expect_equal(c(t2[inds]),rep(0,sum(inds)),tolerance=10^-6,check.attributes=FALSE)
  
  #fit model with diagonal constraint
  modelUncentered=MPLE(X,zUncentered,A,nLambda=101,centered=FALSE,BIC=TRUE,select=TRUE,constraint="diagonal",NR=TRUE)
  
  t1=uncenteredLogLikGrad(cbind(rep(0,14),modelUncentered$MPLE),X,zUncentered,A,constraint="diagonal",thetaInds)
  subModel=modelUncentered$subModel
  
  inds=c(rep(0,14),modelUncentered$subModel)!=0
  t2=uncenteredLogLikGrad(cbind(rep(0,14),modelUncentered$subModel),X,zUncentered,A,constraint="diagonal",thetaInds)
  
  #hopefully near 0
  expect_equal(c(t1),rep(0,length(t1)),tolerance=10^-6,check.attributes=FALSE)
  expect_equal(c(t2[inds]),rep(0,sum(inds)),tolerance=10^-6,check.attributes=FALSE)
})

