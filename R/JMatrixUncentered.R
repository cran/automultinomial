#' Uncentered J matrix
#' 
#' Used for a variance estimate accounting for correlation between sites or for sandwich estimate for independent sites.
#' @param theta a matrix of coefficients
#' @param X a design matrix
#' @param z a matrix of responses
#' @param A a list of adjacency matrices
#' @param constraint "symmetric" or "diagonal"
#' @param thetaInds list of useful indices
#' @param innerIndices the 1-dimensional indices of points within the grid
#' @keywords uncentered, autologistic, variance
#' @return The J matrix for the sandwich variance estimator.
#' @examples
#' #JMatrix=JMatrixUnCentered(theta,X,z,A,innerIndices)
JMatrixUncentered<-function(theta,X,z,A,constraint,thetaInds,innerIndices){
  k=dim(z)[2]
  nBeta=dim(X)[2]
  nEta=k-1
  n=dim(X)[1]
  dldthetaMat=NULL
  #print(dim(theta))
  fullX=as.matrix(cbind(X,A[[1]]%*%z[,2:k]))
  mu=multProb(fullX,theta)
  #print(mu)
  for (i in 1:k){
    dldthetaMat=cbind(dldthetaMat,fullX*(z[innerIndices,i]-mu[,i]))
  }
  
  #get correlated entries
  #newA=A[[1]]
  #newA=Reduce("+",newA)
  #for (i in 1:n){
  #  newA[i,i]=newA[i,i]+1
  #}
  #enforce constraint
  if (constraint=="diagonal" | constraint=="symmetric"){
    transform=toConstraintDerivative(theta,X,z,A,constraint)
    for (i in 1:n){
      temp=NULL
      for (j in 1:k){
        start1=(j-1)*(nBeta+nEta)+nBeta+1
        end1=j*(nEta+nBeta)
        #print(start1)
        #print(end1)
        #print(dldthetaMat[i,])
        #print(dldthetaMat[i,start1:end1])
        temp=c(temp,dldthetaMat[i,start1:end1])
      }
      #print(transform)
      #print(temp)
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
  jMat=Matrix::t(dldthetaMat)%*%A[[1]][,innerIndices]%*%dldthetaMat+Matrix::t(dldthetaMat)%*%dldthetaMat
  return(jMat)
}