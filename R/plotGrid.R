#' Visualize response values on a grid
#'
#' Takes in an nxk matrix or a length n vector containing the response type of each observation
#' and creates an image of the response for a specified grid size
#' @param z nxk matrix or length n vector of responses, where z contains the responses of the grid in column-major order
#' @param m the number of rows in the grid
#' @param n the number of columns in the grid
#' @return void
#' @export
#' @examples
#' #plotGrid(z,m,n)
plotGrid<-function(z,m,n){
  zMat=z
  if (is.vector(z)){
    z=as.numeric(factor(z))
    cats=unique(z)
    zMat=matrix(0,length(z),length(cats))
    for (i in 1:dim(zMat)[2]){
      zMat[i,z[i]]=1
    }
  }
  val=rep(0,m*n)
  for (i in 1:dim(zMat)[2]){
    val=val+i*zMat[,i]
  }
  mat=matrix(val,m,n)
  colors=c("darkolivegreen","darkkhaki","chocolate4","chartreuse4","darkslategrey3","blue","deepskyblue","white","cadetblue2")
  graphics::image(t(mat)[,m:1],col=colors[1:dim(zMat)[2]])
}

