#' Get indices of relevant component of theta matrix
#'
#' Gets indices of beta, eta, diagonal of eta, and off diagonal elements of eta
#' @param theta matrix of coefficients
#' @param X design matrix
#' @param z matrix of multinomial responses
#' @param A list of adjacency matrices
#' @keywords centered
#' @return list of index vectors
#' @export
#' @examples 
#' #
indsFunction<-function(theta,X,z,A){
  #the indices of beta in the (nBeta+nEta)xk theta matrix
  betaInds=NULL
  #the indices of eta in the (nBeta+nEta)xk theta matrix
  etaInds=NULL
  
  #the diagonal elements of the eta matrix in the reference category parameterization
  etaDiag=NULL
  
  #lower and upper off diagonal elements of the eta matrix in reference category parameterization
  etaOffDiagLower=NULL
  etaOffDiagUpper=NULL
  
  k=dim(z)[2]
  nBeta=dim(X)[2]
  nEta=length(A)*(k-1)
  #indices of the reference category etas in the theta matrix
  etaRefInds=rep(FALSE,length(theta))
  etaRefInds[(nBeta+1):(nBeta+nEta)]=TRUE
  theta=matrix(theta,ncol=k)
  betaInds=c(rep(TRUE,nBeta),rep(FALSE,nEta))
  betaInds=rep(betaInds,k)
  etaInds=!betaInds
  etaDiag=etaInds
  etaOffDiagLower=etaInds
  etaOffDiagUpper=etaInds
  for (i in 1:length(A)){
    for (j in 1:k){
      for (m in 2:k){
        index=(nBeta+nEta)*(j-1)+nBeta+(i-1)*(k-1)+(m-1)
        etaDiag[index]=(m==j)
        etaOffDiagLower[index]=(j>1 & m>j)
        etaOffDiagUpper[index]=(j>1 & j>m)
      }
    }
  }
  return(list(betaInds=betaInds,etaInds=etaInds,etaDiag=etaDiag,etaOffDiagLower=etaOffDiagLower,etaOffDiagUpper=etaOffDiagUpper,etaRefInds=etaRefInds))
}


toConstraint<-function(theta,X,z,A,constraint){
  k=dim(z)[2]
  theta=matrix(theta,ncol=k)
  nBeta=dim(X)[2]
  nEta=length(A)*(k-1)
  etaMat=theta[(nBeta+1):(nBeta+nEta),,drop=FALSE]
  etaMat[,1]=0
  if (constraint=="diagonal"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=diag(diag(etaSubMat),nrow=dim(etaSubMat)[2],ncol=dim(etaSubMat)[2])
    }
    theta[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  if (constraint=="symmetric"){
    for (i in 1:length(A)){
      etaSubMat=as.matrix(etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k])
      etaSubMat=as.matrix(Matrix::forceSymmetric(etaSubMat))
      etaMat[((i-1)*(k-1)+1):(i*(k-1)),2:k]=etaSubMat
    }
    theta[(nBeta+1):(nBeta+nEta),]=etaMat
  }
  theta=c(theta)
  return(theta)
}

toConstraintDerivative<-function(theta,X,z,A,constraint){
  k=dim(z)[2]
  theta=matrix(theta,ncol=k)
  nBeta=dim(X)[2]
  nEta=length(A)*(k-1)
  
  transform=matrix(0,(k-1)^2,(k-1)^2)
  count1=1
  inds=matrix(1:((k-1)*(k-1)),k-1,k-1)
  #transformed column
  for (m1 in 1:(k-1)){
    #transformed row
    for (m2 in 1:(k-1)){
      if (m1==m2){
        transform[count1,count1]=1
      }
      if (constraint=="symmetric"){
        if (m1>m2){
          transform[count1,inds[m1,m2]]=1
          transform[count1,inds[m2,m1]]=1
        }
      }
      count1=count1+1
    }
  }
  
  transformTemp=transform
  if (length(A)>1){
    for (i in 1:(length(A)-1)){
      transform=rbind(transform,transformTemp)
    }
  }
  transform=cbind(matrix(0,dim(transform)[1],k-1),transform)
  transform=rbind(matrix(0,k-1,dim(transform)[2]),transform)
  return(transform)
}

#get array index from vector index where the vector corresponds to the strictly above diagonal elements
#of matrix, stored as a vector going down along columns (ie as c(matrix)[aboveDiagonal==TRUE])
getArrayInd<-function(slot){
  col=2
  row=1
  if (slot==1){
    return(c(row,col))
  }
  for (i in 1:(slot-1)){
    if (row+1<col){
      row=row+1
    }
    else{
      row=1
      col=col+1
    }
  }
  return(c(row,col))
}

#given the row and column position of an element on the lower triangle of a matrix, return 
#the slot of that element as it appears in the vector (c(matrix)[belowDiagonal==TRUE])
getVecInd<-function(size,slot){
  ind=1
  row=2
  col=1
  rowEnd=slot[1]
  colEnd=slot[2]
  if (rowEnd==row & colEnd==col){
    return(ind)
  }
  while(rowEnd!=row | colEnd!=col){
    ind=ind+1
    if (row+1<=size){
      row=row+1
    }
    else{
      col=col+1
      row=col+1
    }
  }
  return(ind)
}
slotUpper_to_slotLower<-function(size,slot){val=getArrayInd(slot);return(getVecInd(size,c(val[2],val[1])))}

#we have the upper triangular variance components already in the matrix, now copy these over for the lower triangular part
varianceTransformation<-function(theta,covariance){
  indsUpper=rep(FALSE,length(theta))
  indsLower=rep(FALSE,length(theta))
  k=dim(theta)[2]+1
  nBeta=dim(theta)[1]-(k-1)
  #get upper/lower triangular parts
  #cycle through theta
  count=1
  for (i in 1:dim(theta)[2]){
    for (j in 1:dim(theta)[1]){
      etaRow=j-nBeta
      etaCol=i
      #if we're in the eta section
      if (etaRow>0){
        #lower triangle
        if (etaRow>etaCol){
          indsLower[count]=TRUE
        }
        if (etaRow<etaCol){
          indsUpper[count]=TRUE
        }
      }
      count=count+1
    }
  }
  
  #now we have a logical vector, but we really want the actual indices
  indsUpper=which(indsUpper)
  indsLower=which(indsLower)
  
  #now we can update the covariance matrix
  for (i in 1:length(indsUpper)){
    upperIndex=indsUpper[i]
    lowerIndex=indsLower[slotUpper_to_slotLower(k-1,i)]
    covariance[lowerIndex,]=covariance[upperIndex,]
    covariance[,lowerIndex]=covariance[,upperIndex]
  }
  
  for (i in 1:length(indsUpper)){
    #correlation of 1 b/w upper and lower elements due to symmetry
    upperIndex=indsUpper[i]
    lowerIndex=indsLower[slotUpper_to_slotLower(k-1,i)]
    covariance[lowerIndex,upperIndex]=covariance[upperIndex,upperIndex]
    covariance[upperIndex,lowerIndex]=covariance[upperIndex,upperIndex]
  }
  return(covariance)
}

boundaryIndices<-function(nRows,nCols){
  
  #get site indices of boundary sites
  #don't count the corners twice
  nBoundary=nRows+nRows+(nCols-2)+(nCols-2)
  boundVec=1:nBoundary
  
  #the first column of the grid is made up of boundary points
  boundVec[1:nRows]=1:nRows
  
  column=2
  count=nRows
  #get top and bottom points
  for (i in 2:(nCols-1)){
    count=count+1
    boundVec[count]=(column-1)*nRows+1
    count=count+1
    boundVec[count]=column*nRows
    column=column+1
  }
  count=count+1
  boundVec[count:nBoundary]=nRows*(nCols-1)+1:nRows
  return(boundVec)
}