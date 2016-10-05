#' Spatial adjacency matrix
#'
#' Takes as input a matrix gridMat with the following format: 0 means a site is included in the study area, and nonzero values mean the site is not included.
#' Returns an adjacency matrix for the study sites, where sites at a Manhattan distance of exactly kNN are considered neighbors. 
#' @param m height of the grid
#' @param n width of the grid
#' @param boundary surround grid by an external boundary, defaults to FALSE
#' @keywords adjacency, spatial, neighborhood, grid
#' @export
#' @return a list containing a sparse adjacency matrix for a square lattice, and the vector of the points not on the edge
#' @examples 
#' mat=adjMat(10,10)
adjMat<-function(m,n,boundary=FALSE){
  nRows=m
  nCols=n
  if (boundary){
    nRows=nRows+2
    nCols=nCols+2
  }
  
  n1=nRows*nCols
  A=Matrix::Matrix(0,n1,n1)
  
  g<-igraph::make_lattice(c(nRows,nCols))
  A=igraph::get.adjacency(g)
  
  #if we have boundary conditions, remove rows corresponding to entries past the edge of the grid
  all=1:n1
  innerIndices=all
  if (boundary){
    innerIndices=which(!(all%in%boundaryIndices(nRows,nCols)))
  }
  return(list(A=A,innerIndices=innerIndices))
}