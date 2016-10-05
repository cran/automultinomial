#' J matrix for centered variance estimate
#'
#' Computes J matrix for the centered variance estimate
#' @param theta coefficients
#' @param X design matrix 
#' @param z response matrix
#' @param A list of adjacency matrices
#' @param constraint "symmetric" or "diagonal"
#' @param thetaInds list of relevant indices in theta: betaInds, etaDiag, ...
#' @keywords adjacency, spatial, neighborhood, grid
#' @return A sparse adjacency matrix
#' @examples 
#' #JMatrix=JMatrixCentered(theta,X,z,A)
JMatrixCentered<-function(theta,X,z,A,constraint,thetaInds){
  k=dim(z)[2]
  n=dim(z)[1]
  nEta=length(A)*(k-1)
  nBeta=dim(X)[2]
  theta=matrix(theta,nrow=nBeta+nEta)
  newX=totalX(theta,X,z,A)
  dldthetaMat=NULL
  for (i in 1:k){
    dldthetaMat=cbind(dldthetaMat,newX*z[,i])
  }
  
  
  eta=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  
  mu=multProb(X,theta[1:nBeta,,drop=FALSE])
  
  
  
  count=0
  gMat=matrix(0,n,k*nBeta)
  gradient=matrix(0,n,k*nBeta)
  #go through each coefficient
  #outer: neighborhood of column in X matrix
  #get the part of the gradient from the trace
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
      for (i in 1:k){
        #get the eta corresponding to A[[m]] and mu[,k] for response i
        etaCurrent=eta[count,i]
        inds=z[,i]==1
        gradient=gradient+etaCurrent*inds*gMat
      }
    }
  }
  
  temp=NULL
  for (i in 1:k){
    start=(i-1)*nBeta+1
    end=i*nBeta
    temp=cbind(temp,gradient[,start:end],replicate(nEta,rep(0,n)))
  }
  dldthetaMat=dldthetaMat+temp
  
  #get the part of the gradient from the exp sum
  gradBeta=matrix(0,n,k*nBeta)
  gradEta=matrix(0,n,(k)*nEta)
  denominator=apply(exp(newX%*%theta),1,sum)
  for (j in 1:k){
    start=(j-1)*nBeta+1
    end=j*nBeta
    #gradient of exp(newX%*%theta_m) wrt beta_j
    etaStart=(j-1)*nEta+1
    etaEnd=j*nEta
    for (m in 1:k){
      tempExp=as.vector(exp(newX%*%theta[,m]))
      if (m==j){
        val=X*tempExp/denominator
        val=as.matrix(val)
        gradBeta[,start:end]=gradBeta[,start:end]+val
        
        val=newX[,(nBeta+1):(nBeta+nEta)]*tempExp/denominator
        val=as.matrix(val)
        gradEta[,etaStart:etaEnd]=gradEta[,etaStart:etaEnd]+val
      }
      
      #now go through each neighborhood/eta
      second=rep(0,nBeta)
      count=0
      for (q in 1:length(A)){
        for (r in 2:k){
          count=count+1
          tempW=mu[,r]*((r==j)-mu[,j])
          etaTemp=eta[count,m]
          val=-(A[[q]]%*%(tempW*X))*1/denominator*etaTemp*tempExp
          val=as.matrix(val)
          gradBeta[,start:end]=gradBeta[,start:end]+val
        }
      }
    }
  }
  
  temp=NULL

  
  for (i in 1:k){
    start1=(i-1)*nBeta+1
    end1=i*nBeta
    start2=(i-1)*nEta+1
    end2=i*nEta
    temp=cbind(temp,gradBeta[,start1:end1],gradEta[,start2:end2])
  }
  dldthetaMat=dldthetaMat+-temp
  
  
  #enforce constraint
  if (constraint=="diagonal" | constraint=="symmetric"){
    transform=toConstraintDerivative(theta,X,z,A,constraint)
    for (i in 1:n){
      temp=NULL
      for (j in 1:k){
        start1=(j-1)*(nBeta+nEta)+nBeta+1
        end1=j*(nEta+nBeta)
        temp=c(temp,dldthetaMat[i,start1:end1])
      }
      temp=transform%*%temp
      for (j in 1:k){
        start1=(j-1)*(nBeta+nEta)+nBeta+1
        end1=j*(nEta+nBeta)
        start2=(j-1)*nEta+1
        end2=j*nEta
        dldthetaMat[i,start1:end1]=temp[start2:end2]
      }
    }
  }
  
  #get TRUE/FALSE values of neighborhood structure for each obs.
  newA=A
  for (i in 1:length(A)){
    newA[[i]]=(newA[[i]]!=0)
  }
  newA=Reduce("+",newA)
  newA=(newA!=FALSE)*1
  newA=newA
  jMat=Matrix::t(dldthetaMat)%*%newA%*%dldthetaMat+Matrix::t(dldthetaMat)%*%dldthetaMat
  return(jMat)
}