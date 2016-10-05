traceHessian<-function(theta,X,z,A){
  k=dim(z)[2]
  n=dim(z)[1]
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  theta=matrix(theta,nrow=nBeta+nEta)
  #newX=totalX(theta,X,z,A)
  
  eta=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])
  
  
  
  count=0
  hessian=matrix(0,k*nBeta,k*nBeta)
  etaHess=matrix(0,k*nBeta,k*nEta)
  for (j in 1:(length(A))){
    
    #the category of this X1 matrix column
    for (m in 2:k){
      count=count+1
      #gradient of likelihood wrt beta_l
      gMat=matrix(0,n,k*nBeta)
      for (l in 1:k){
        start1=(l-1)*nBeta+1
        end1=l*nBeta
        
        #jacobian wrt beta_s of gradient wrt beta_r
        for (r in l:k){
          start2=(r-1)*nBeta+1
          end2=r*nBeta
          
          weights=mu[,m]*((m==r)-mu[,r])*((m==l)-mu[,l])
          weights=weights+mu[,m]*(-mu[,l]*((l==r)-mu[,r]))
          weights=-weights
          subHess=matrix(0,nBeta,nBeta)
          for (s in 1:k){
            inds=(z[1:n,s]==1)
            etaTemp=eta[count,s]
            Aweights=apply(A[[j]][inds,],2,sum)
            Aweights=Aweights*weights
            Aweights=Aweights*etaTemp
            subHess=subHess+Matrix::crossprod(X,Aweights*X)
          }
          hessian[start1:end1,start2:end2]=hessian[start1:end1,start2:end2]+subHess
        }
        
        #jacobian wrt eta of gradient wrt beta
        if (l==m){
          gMat[,start1:end1]=-as.matrix(A[[j]]%*%(mu[,l]*(1-mu[,l])*X))
        }
        if (l!=m){
          gMat[,start1:end1]=-as.matrix(A[[j]]%*%((mu[,l]*-mu[,m])*X))
        }
      }
      
      
      indStart=count
      inds=c(indStart,indStart+(1:(k-1))*nEta)
      #print(inds)
      etaHess[,inds]=etaHess[,inds]+t(gMat)%*%z
    }
  }
  hessian=as.matrix(Matrix::forceSymmetric(hessian))
  p=nBeta+nEta
  full=matrix(0,k*p,k*p)
  
  indStart=1:nBeta
  betaInds=indStart
  for (i in 1:(k-1)){
    betaInds=c(betaInds,indStart+i*p)
    #print(inds)
  }
  full[betaInds,betaInds]=hessian
  
  indStart=(nBeta+1):p
  etaInds=indStart
  for (i in 1:(k-1)){
    etaInds=c(etaInds,indStart+i*p)
  }
  #print(inds)
  #print(etaInds)
  #print(betaInds)
  full[betaInds,etaInds]=etaHess
  full[etaInds,betaInds]=t(etaHess)
  return(full)
}


expSumHessian<-function(theta,X,z,A){
  
  k=dim(z)[2]
  n=dim(z)[1]
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  theta=matrix(theta,nrow=nBeta+nEta)
  newX=totalX(theta,X,z,A)
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])
  muTot=multProb(newX,theta)
  qBeta=rep(0,nBeta*k)
  gradEta=rep(0,nEta*k)
  eta=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  p=nBeta*k
  pTot=(nBeta+nEta)*k
  betaHess=matrix(0,pTot,pTot)
  
  denominator=apply(exp(newX%*%theta),1,sum)
  tempExp=exp(newX%*%theta)
  #suppose we took a derivative wrt beta_j
  for (j in 1:k){
    start1=(j-1)*(nBeta+nEta)+1
    end1=(j-1)*(nBeta+nEta)+nBeta
    
    #now we'll do a Jacobian with respect to beta_i
    for (i in 1:k){
      start2=(i-1)*(nBeta+nEta)+1
      end2=(i-1)*(nBeta+nEta)+nBeta
      
      #get the Qi term
      qBetai=rep(0,nBeta)
      for (m in 1:k){
        #print("Hello")
        #tempExpM=as.vector(exp(newX%*%theta[,m]))
        tempExpM=tempExp[,m]
        qBetai=qBetai+(m==i)*X*tempExpM
        
        count=0
        for (r in 1:length(A)){
          for (s in 2:k){
            count=count+1
            weight1=mu[,s]*((i==s)-mu[,i])
            etaTemp=eta[count,m]
            qBetai=qBetai+(-etaTemp*A[[r]]%*%(weight1*X))*tempExpM
          }
        }
      }
      
      #now the qi part is done
      
      
      #the derivative wrt beta_j was a sum, so we have to differentiate each term
      for (l in 1:k){
        #get first term
        #print("Hello")
        #tempExpL=as.vector(exp(newX%*%theta[,l]))
        tempExpL=tempExp[,l]
        t1=-tempExpL/(denominator^2)*qBetai
        temp=rep(0,nBeta)
        weight=tempExpL/denominator
        temp=temp+(i==l)*X*weight
        
        #get f function
        count=0
        for (r in 1:length(A)){
          for (s in 2:k){
            count=count+1
            weight1=mu[,s]*((i==s)-mu[,i])
            etaTemp=eta[count,l]
            temp=temp+(-etaTemp*A[[r]]%*%(weight1*X))*weight
          }
        }
        #weighted gradient of mu_l wrt beta_i
        t1=t1+temp
        
        #get d2/dbetadeta
        if (l==j){
          startEta=(j-1)*(nBeta+nEta)+nBeta+1
          endEta=(j-1)*(nBeta+nEta)+nBeta+nEta
          first=t1
          second=newX[,(nBeta+1):(nBeta+nEta)]
          val=as.matrix(Matrix::crossprod(first,second))
          betaHess[start2:end2,startEta:endEta]=betaHess[start2:end2,startEta:endEta]+val
          betaHess[startEta:endEta,start2:end2]=t(betaHess[start2:end2,startEta:endEta])
          
          weight=muTot[,j]
          count=0
          for (r in 1:length(A)){
            for (s in 2:k){
              count=count+1
              ind=startEta+count-1
              weight1=mu[,s]*((i==s)-mu[,i])
              val=(-A[[r]]%*%(weight1*X))*weight
              val=apply(val,2,sum)
              betaHess[ind,start2:end2]=betaHess[ind,start2:end2]+val
              betaHess[start2:end2,ind]=betaHess[ind,start2:end2]
            }
          }
        }
        #get second term
        t2=(j==l)*X
        count=0
        for (r in 1:length(A)){
          for (s in 2:k){
            count=count+1
            weight1=mu[,s]*((j==s)-mu[,j])
            etaTemp=eta[count,l]
            t2=t2+(-etaTemp*A[[r]]%*%(weight1*X))
          }
        }
        val=Matrix::crossprod(t2,t1)
        betaHess[start1:end1,start2:end2]=betaHess[start1:end1,start2:end2]+as.matrix(val)
      }
    }
  }
  #second term
  for (i in 1:k){
    start1=(i-1)*(nBeta+nEta)+1
    end1=(i-1)*(nBeta+nEta)+nBeta
    for (j in 1:k){
      start2=(j-1)*(nBeta+nEta)+1
      end2=(j-1)*(nBeta+nEta)+nBeta
      test1=matrix(0,nBeta,nBeta)
      test2=test1
      for (l in 1:k){
        weight=muTot[,l]
        
        count=0
        for (r in 1:length(A)){
          for (s in 2:k){
            count=count+1
            inds=apply(weight*A[[r]],2,sum)
            etaTemp=eta[count,l]
            mu_s=mu[,s]
            mu_i=mu[,i]
            mu_j=mu[,j]
            
            weights=mu_s*((i==s)-mu_i)*((j==s)-mu_j)
            weights=weights+mu_s*(-mu_j*((i==j)-mu_i))
            Aweights=as.vector(inds*weights)
            
            subHess=-Matrix::crossprod(X,Aweights*X)*etaTemp
            betaHess[start1:end1,start2:end2]=betaHess[start1:end1,start2:end2]+subHess
          }
        }
      }
    }
  }
  for (i in 1:k){      
    start1=(i-1)*(nBeta+nEta)+nBeta+1
    end1=(i-1)*(nBeta+nEta)+nBeta+nEta
    for (j in 1:k){
      start2=(j-1)*(nBeta+nEta)+nBeta+1
      end2=(j-1)*(nBeta+nEta)+nBeta+nEta
      first=newX[,(nBeta+1):(nBeta+nEta)]
      weight=muTot[,i]*((i==j)-muTot[,j])
      val=Matrix::crossprod(first,weight*first)
      betaHess[start1:end1,start2:end2]=betaHess[start1:end1,start2:end2]+val
    }
  }
  betaHess=as.matrix(Matrix::forceSymmetric(betaHess))
  return(betaHess)
}
#' Centered model hessian
#' 
#' Computes a hessian for the centered model, taking coefficients, a design matrix, a response matrix, and a list of adjacency matrices 
#' @param theta a matrix of coefficients with dimension (p+(k-1)*length(A)) by k; frequently A will have length 1
#' @param X nxp design matrix
#' @param z nxk response matrix
#' @param A a list of adjacency matrices
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @keywords centered, gradient
#' @return a gradient vector
#' @export
#' @examples 
#' #g=centeredLogLikHess(theta,X,z,A)
centeredLogLikHess<-function(theta,X,z,A,constraint,thetaInds){
  if (constraint != "none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  val1=traceHessian(theta,X,z,A)
  val2=expSumHessian(theta,X,z,A)
  val=val1-val2
  if (constraint=="none"){
    return(val)
  }
  if (constraint=="symmetric" | constraint=="diagonal"){
    etaHess=val[thetaInds$etaInds,thetaInds$etaInds,drop=FALSE]
    transform=toConstraintDerivative(theta,X,z,A,constraint)
    etaHess=transform%*%etaHess%*%t(transform)
    betaEtaHess=val[thetaInds$etaInds,thetaInds$betaInds,drop=FALSE]
    betaEtaHess=transform%*%betaEtaHess
    val[thetaInds$etaInds,thetaInds$etaInds]=etaHess
    val[thetaInds$etaInds,thetaInds$betaInds]=betaEtaHess
    val[thetaInds$betaInds,thetaInds$etaInds]=t(betaEtaHess)
  }
  return(as.matrix(Matrix::forceSymmetric(val)))
}