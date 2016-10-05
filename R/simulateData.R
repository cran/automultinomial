#' Simulate data
#' 
#' Simulates data from the centered or uncentered binomial or multinomial distribution
#' @param theta a vector of coefficients for binary response, or a matrix of coefficients for multicategory response
#' @param X n by p design matrix
#' @param m width of the grid
#' @param n height of the grid
#' @param boundary (uncentered only) -1 for no boundary conditions, i>=0 to surround grid with class i observations 
#' @param centered logical, generates data from centered model if TRUE and uncentered if FALSE
#' @param k number of response categories
#' @param burn number of burn-in iterations for the Gibbs sampler
#' @param draws number of simulated samples to return
#' @keywords simulate
#' @export
#' @return an nxk matrix of responses for draws=1, and a length draws list of nxk response matrices for draws>1
#' @examples
#' 
#' #generate outcomes on a grid using example coefficient values and design matrix
#' set.seed(42)
#' #for a multinomial response, beta and eta will have more than one column
#' #for a binary response, beta and eta will only have one column
#' #here, we the response variable takes 3 possible values, so beta and eta have two columns
#' beta=cbind(c(0,1.25,-1,2,-0.8),c(-3,0.25,-5,1,-0.4))
#' 
#' #setting the eta coefficients
#' eta=cbind(c(0.6,0),c(0,0.6))
#' 
#' #these are the final coefficient values we'll use to generate the data
#' theta=rbind(beta,eta)
#' 
#' #X matrix with 5 predictors and 900 observations for a 10x10 grid
#' X=cbind(rep(1,900),replicate(4,rnorm(900)))
#' 
#' #generate data for 30x30 grid
#' z=simulateData(theta,X,30,30,centered=TRUE,k=3,burn=1)
simulateData<-function(theta,X,m,n,boundary=-1,centered=FALSE,k=2,burn=40,draws=1){
  boundFlag=(boundary>-1 & !centered)
  #adjacency matrix for full grid
  A=NULL
  innerMat=matrix(1:(m*n),m,n)
  if (boundFlag){
    m=m+2
    n=n+2
    innerMat=rbind(rep(0,n-2),innerMat,rep(0,n-2))
    innerMat=cbind(rep(0,m),innerMat,rep(0,m))
  }
  
  A=adjMat(m,n,FALSE)$A
  
  rmultinom2<-function(probs){
    zOutcomes=matrix(0,dim(probs)[1],dim(probs)[2])
    outcomes=Hmisc::rMultinom(probs,1)
    for (i in 1:dim(probs)[1]){
      zOutcomes[i,outcomes[i,1]]=1
    }
    return(zOutcomes)
  }
  
  thin=1
  nObs=m*n
  count=0
  nEta=(k-1)
  nBeta=dim(as.matrix(theta))[1]-nEta
  
  
  #boundary and nonBoundary
  all=NULL
  boundaryPoints=NULL
  nonBoundaryPoints=NULL
  if (boundFlag){
    all=1:nObs
    boundaryPoints=boundaryIndices(m,n)
    nonBoundaryPoints=which(!(all%in%boundaryPoints))
  }

  
  #put theta in mean 0 parameterization
  theta=cbind(rep(0,nEta+nBeta),theta)
  theta=theta-apply(theta,1,mean)
  
  mu=multProb(X,theta[1:nBeta,])
  z=matrix(0,nObs,k)
  count=0
  
  if (boundFlag){
    for (j in nonBoundaryPoints){
      count=count+1
      z[j,]=stats::rmultinom(1,1,mu[count,])
    }
    z[boundaryPoints,boundary+1]=1
  }
  if (!boundFlag){
    for (j in nObs){
      z[j,]=stats::rmultinom(1,1,mu[j,])
    }
  }
  
  Agraph=igraph::graph_from_adjacency_matrix(A)
  
  #firstVertices and secondVertices are conditionally independent
  firstVertices=igraph::bipartite.mapping(Agraph)$type
  secondVertices=!firstVertices
  firstVertices=which(firstVertices)
  secondVertices=which(secondVertices)
  
  #get rid of boundary rows
  firstVerticesSubset=firstVertices
  secondVerticesSubset=secondVertices
  if (boundFlag){
    firstVerticesSubset=firstVertices[!(firstVertices%in%boundaryPoints)]
    secondVerticesSubset=secondVertices[!(secondVertices%in%boundaryPoints)]
  }
  
  A_vertices=list(firstVertices,secondVertices)
  A_verticesSubset=list(firstVerticesSubset,secondVerticesSubset)
  X_vertices=list(innerMat[firstVerticesSubset],innerMat[secondVerticesSubset])
  
  XjList=list(X[X_vertices[[1]],,drop=FALSE],X[X_vertices[[2]],,drop=FALSE])
  Alist=list(A[A_verticesSubset[[1]],A_vertices[[2]]],A[A_verticesSubset[[2]],A_vertices[[1]]])
  
  for (i in 1:(burn)){
    if (centered){
      for (j in 1:2){
        jOther=2-j+1        
          
        Xj=cbind(XjList[[j]],Alist[[j]]%*%(z[A_vertices[[jOther]],2:k]-mu[A_vertices[[jOther]],2:k]))
        probs=multProb(Xj,theta)
        z[A_verticesSubset[[j]],]=rmultinom2(probs)

#         for (m in 1:length(vertices[[j]])){
#           z[vertices[[j]][m],]=rmultinom(1,1,probs[m,])
#         }
      }
    }
    if (!centered){
      for (j in 1:2){
        jOther=2-j+1
        Xj=cbind(XjList[[j]],Alist[[j]]%*%(z[A_vertices[[jOther]],2:k]))
        probs=multProb(Xj,theta)
        z[A_verticesSubset[[j]],]=rmultinom2(probs)
        
#        for (m in 1:length(vertices[[j]])){
#           z[vertices[[j]][m],]=rmultinom(1,1,probs[m,])
#         }
      }
    }
  }
  zList=list()
  if (draws==1){return(z)}
  zList[[1]]=z
  for (h in 2:draws){
    for (i in 1:thin){
      if (centered){
        for (j in 1:2){
          jOther=2-j+1
          
          Xj=cbind(XjList[[j]],Alist[[j]]%*%(z[A_vertices[[jOther]],2:k]-mu[A_vertices[[jOther]],2:k]))
          
          probs=multProb(Xj,theta)
          z[A_verticesSubset[[j]],]=rmultinom2(probs)
          
#           for (m in 1:length(vertices[[j]])){
#             z[vertices[[j]][m],]=rmultinom(1,1,probs[m,])
#           }
        }
      }
      if (!centered){
        for (j in 1:2){
          jOther=2-j+1
          
          Xj=cbind(XjList[[j]],Alist[[j]]%*%(z[A_vertices[[jOther]],2:k]))
          
          probs=multProb(Xj,theta)
          z[A_verticesSubset[[j]],]=rmultinom2(probs)
          
#           for (m in 1:length(vertices[[j]])){
#             z[vertices[[j]][m],]=rmultinom(1,1,probs[m,])
#           }
        }
      }
    }
    zList[[h]]=z
  }
  return(zList)
}