H_stat_2s_2d <- function(x,y,tol=1e-6){
  if(is.list(x)){
    if(length(x) != 2) stop("The sample has to be two dimensional")
    x <-  matrix(unlist(x),ncol = 2)
  } 
  
  if(is.list(y)){
    if(length(y) != 2) stop("The sample has to be two dimensional")
    y <-  matrix(unlist(y),ncol = 2)
  }
  
  # Compute Boundary Values M
  M <- round(max(max(x),max(y)))+2.5
  
  # compute the omnidirectional of the empirical function of x
  x1 <- sort(x[,1])
  x2 <- sort(x[,2])
  n1 <- length(x1)
  n2 <- length(x2)
  # also include the boundary vertices (step 2)
  x_rectl <- list(c(x1-tol,x1,M),
                  c(x2-tol,x2,M))
  z_x <- multiecdf(x,z=x_rectl)
  n1_rectl <- length(x_rectl[[1]])
  z_x <- matrix(z_x, ncol = n1_rectl, byrow = T)
  x_jumps <- NULL
  x_jumps <- NULL
  
  x_rectl2 <- list(sort(c(x1,M)),
                  sort(c(x2,M)))
  z_yx <- multiecdf(y,z=x_rectl2)
  z_yx <- matrix(z_yx, ncol = n1+1, byrow = T)
  for (xi in 1:n1) {
    for (yi in 1:n2) {
      if( (z_x[2*yi,2*xi]!=z_x[2*yi,2*xi-1])&&(z_x[2*yi,2*xi]!=z_x[2*yi-1,2*xi]) ) {
        z_1_jump <- z_x[2*yi,2*xi]
        z_2_jump <- z_x[2*yi-1,2*xi]
        z_3_jump <- z_x[2*yi,2*xi-1]
        z_4_jump <- z_x[2*yi-1,2*xi-1]
        z_y <- z_yx[yi,xi]
        x_jump <- c(x_rectl[[1]][2*xi],x_rectl[[2]][2*yi])
        # compute the index,location of jumps, projection of vertices
        x_jumps <- rbind(x_jumps,c(x_jump, z_1_jump, z_2_jump, z_3_jump, z_4_jump,z_y))
      }   
    }
    z_1_jump <- z_x[2*n2+1,2*xi]
    z_4_jump <- z_x[2*n2+1,2*xi-1]
    if( z_1_jump!=z_4_jump ){
      z_y <- z_yx[n2+1,xi]
      x_jump <- c(x_rectl[[1]][2*xi],M)
      x_jumps <- rbind(x_jumps,c(x_jump, z_1_jump, rep((z_1_jump+z_4_jump)/2,2), z_4_jump,z_y))
    }
  }
  for(yi in 1:n2) {
    z_1_jump <-  z_x[2*yi,2*n1+1]
    z_4_jump <-  z_x[2*yi-1,2*n1+1]
    if( z_1_jump!=z_4_jump ) {
      z_y <- z_yx[yi,n1+1]
      x_jump <- c(M,x_rectl[[2]][2*yi])
      x_jumps <- rbind(x_jumps,c(x_jump, z_1_jump, rep((z_1_jump+z_4_jump)/2,2),z_4_jump,z_y))
    }
  }
  projection_x <- NULL
  # no tie
  for (i in 1:nrow(x_jumps)) {
    if(x_jumps[i,3]>x_jumps[i,7]) projection_x <- rbind(projection_x, c(x_jumps[i,1:2]+x_jumps[i,3],x_jumps[i,1:3]))
    if((x_jumps[i,7]>x_jumps[i,6]+1e-14)&&(abs(x_jumps[i,6]-x_jumps[5])>1e-14)) projection_x <- rbind(projection_x, c(x_jumps[i,1:2]+x_jumps[i,6],x_jumps[i,1:2],x_jumps[i,6]))
  }
  
  # compute the omnidirectional of the empirical function of y
  y1 <- sort(y[,1])
  y2 <- sort(y[,2])
  n1 <- length(y1)
  n2 <- length(y2)
  y_rectl <- list(c(y1-tol,y1),
                  c(y2-tol,y2))
  z_y <- multiecdf(y,z=y_rectl)
  n1_rectl <- length(y_rectl[[1]])
  z_y <- matrix(z_y, ncol = n1_rectl, byrow = T)
  
  projection_y <- NULL
  for (xi in 1:n1) {
    for (yi in 1:n2) {
      if( (z_y[2*yi,2*xi]!=z_y[2*yi,2*xi-1])&&(z_y[2*yi,2*xi]!=z_y[2*yi-1,2*xi]) ) {

        y_jump <- c(y_rectl[[1]][2*xi],y_rectl[[2]][2*yi])
        z_1_jump <- z_y[2*yi,2*xi]
        z_2_jump <- z_y[2*yi-1,2*xi]
        z_3_jump <- z_y[2*yi,2*xi-1]
        z_4_jump <- z_y[2*yi-1,2*xi-1]
        projection_y <- rbind(projection_y, c(y_jump+z_2_jump,y_jump, z_2_jump,2))
        if(z_4_jump != z_2_jump){
          projection_y <- rbind(projection_y, c(y_jump+z_4_jump,y_jump, z_4_jump,4))
          if(z_2_jump != z_3_jump){
            projection_y <- rbind(projection_y, c(y_jump+z_3_jump,y_jump, z_3_jump,3))
          }
        }
        projection_y <- rbind(projection_y, c(y_jump+z_1_jump,y_jump, z_1_jump,1))
      }
    }
  }

  
  
  
  ##### compute H
  return(hsearch_Rcpp(projection_x,projection_y))
}
    

