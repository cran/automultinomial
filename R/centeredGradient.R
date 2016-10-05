traceGrad<-function(theta,X,z,A){
  k=dim(z)[2]
  n=dim(z)[1]
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  theta=matrix(theta,nrow=nBeta+nEta)
  newX=totalX(theta,X,z,A)
  first=t(newX)%*%z
  first=c(first)
  
  eta=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  z_etat=z%*%t(eta)
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])
  
  
  count=0
  gMat=matrix(0,n,k*nBeta)
  gradient=rep(0,nBeta*k)
  
  #go through each coefficient
  #outer: neighborhood of column in X matrix
  for (j in 1:(length(A))){
    
    #inner: category of column in X matrix
    #gMat will be the matrix where each row
    #is the gradient wrt beta_l of (A[[j]]*-mu[,m])<-- this is a vector
    for (m in 2:k){
      count=count+1
      
      #category of beta coefficient
      for(l in 1:k){
        start=(l-1)*nBeta+1
        end=l*nBeta
        if (l==m){
          gMat[,start:end]=-as.matrix(A[[j]]%*%(mu[,l]*(1-mu[,l])*X))
        }
        if (l!=m){
          gMat[,start:end]=-as.matrix(A[[j]]%*%((mu[,l]*-mu[,m])*X))
        }
      }
      gradient=gradient+as.vector(t(gMat)%*%z_etat[,count,drop=FALSE])
    }
  }
  gradient=matrix(gradient,ncol=k)
  gradient=rbind(gradient,replicate(k,rep(0,nEta)))
  gradient=first+c(gradient)
  return(gradient)
}

expSumGradient<-function(theta,X,z,A){
  k=dim(z)[2]
  n=dim(z)[1]
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  theta=matrix(theta,nrow=nBeta+nEta)
  newX=totalX(theta,X,z,A)
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])
  gradBeta=rep(0,nBeta*k)
  gradEta=rep(0,nEta*k)
  eta=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  
  #go through each category, take gradient wrt beta_j
  denominator=apply(exp(newX%*%theta),1,sum)
  tempExpMat=exp(newX%*%theta)
  for (j in 1:k){
    start=(j-1)*nBeta+1
    end=j*nBeta
    #gradient of exp(newX%*%theta_m) wrt beta_j
    etaStart=(j-1)*nEta+1
    etaEnd=j*nEta
    for (m in 1:k){
      tempExp=tempExpMat[,m]
      if (m==j){
        val=X*(tempExp/denominator)
        val=apply(val,2,sum)
        gradBeta[start:end]=gradBeta[start:end]+val
        
        val=newX[,(nBeta+1):(nBeta+nEta)]*(tempExp/denominator)
        if (nEta>1){
          val=apply(val,2,sum)
        }
        if (nEta==1){
          val=sum(val)
        }
        gradEta[etaStart:etaEnd]=gradEta[etaStart:etaEnd]+val
      }
      
      #now go through each neighborhood/eta
      second=rep(0,nBeta)
      count=0
      for (q in 1:length(A)){
        for (r in 2:k){
          count=count+1
          tempW=mu[,r]*((r==j)-mu[,j])
          etaTemp=eta[count,m]
          val=-(A[[q]]%*%(tempW*X))*(1/denominator*etaTemp*tempExp)
          val=apply(val,2,sum)
          gradBeta[start:end]=gradBeta[start:end]+val
        }
      }
    }
  }
  
  gradBeta=matrix(gradBeta,ncol=k)
  gradEta=matrix(gradEta,ncol=k)
  gradient=c(rbind(gradBeta,gradEta))
  return(gradient)
}

#' Centered model gradient
#' 
#' Computes a gradient for the centered model, taking coefficients, a design matrix, a response matrix, and a list of adjacency matrices 
#' @param theta a matrix of coefficients with dimension (p+(k-1)*length(A)) by k; frequently A will have length 1
#' @param X nxp design matrix
#' @param z nxk response matrix
#' @param A a list of adjacency matrices
#' @param constraint constraint on eta
#' @param thetaInds relevant indices
#' @param innerIndices only here to match argument list with uncentered case
#' @keywords centered, gradient
#' @export
#' @return a gradient vector
#' @examples 
#' #g=centeredLogLikGrad(theta,X,z,A)
centeredLogLikGrad<-function(theta,X,z,A,constraint,thetaInds,innerIndices=NULL){
  if (constraint!="none"){
    theta=toConstraint(theta,X,z,A,constraint)
  }
  val=traceGrad(theta,X,z,A)-expSumGradient(theta,X,z,A)
  if (constraint=="none"){
    return(val)
  }
  k=dim(z)[2]
  val=matrix(val,ncol=k)
  nBeta=dim(X)[2]
  nEta=length(A)*(k-1)
  etaMat=val[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  etaMat[,1]=0
  if (constraint=="diagonal"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=diag(diag(etaSubMat),nrow=dim(etaSubMat)[2],ncol=dim(etaSubMat)[2])
    }
    val[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  if (constraint=="symmetric"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaSubMat[upper.tri(etaSubMat)]=etaSubMat[upper.tri(etaSubMat)]+t(etaSubMat)[upper.tri(t(etaSubMat))]
      etaSubMat[lower.tri(etaSubMat)]=0
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=etaSubMat
    }
    val[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  return(val)
}