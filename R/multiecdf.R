multiecdf <- function(x,z = NULL){

  if(is.list(x)){
    x <- matrix(unlist(x), byrow=F, ncol=length(x) )
  } else if(is.data.frame(x)){
    x <- as.matrix(x)
  }

  if(is.null(z)){
    return(fastCDFOnSample_Rcpp(p_x_r = t(x)))
  } else {
    if(length(z)!=ncol(x)) stop("The evaluation grid has to have the same dimension as sample")
    return(fastCDF_Rcpp(p_x_r = t(x), p_z_r = z))
  }
}


