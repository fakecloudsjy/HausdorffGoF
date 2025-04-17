multiecdf <- function(x,z = NULL){
  if(is.list(x)){
    x <- matrix(unlist(x), byrow=TRUE, nrow=length(x) )
  }
  y <- rep(1,ncol(x))
  if(is.null(z)){
    return(fastCDFOnSample_wrapper(p_x_r = x, p_y_r = y))
  } else {
    if(length(z)!=nrow(p_x)) stop("The evaluation grid has to have the same dimension as sample")
    return(fastCDF_wrapper(p_x_r = x, p_z_r = z, p_y_r = y))
  }
}