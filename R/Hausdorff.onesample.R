H_test_1s_1d = function(x, CDF, CDFinverse = NULL, pdf = NULL, tol = 1e-10){
  q = H_stat_1s_1d(x, CDF, pdf, tol)
  n = length(x)
  f_b = rep(0,n)
  f_a = rep(0,n)
  if(is.null(CDFinverse)){
    f <- function(lambda){
      function(x){
        CDF(x) - lambda
      }
    }
    if(is.null(pdf)){
      for (i in 1:n) {
        if(i/n-q<0){
          f_a[i] <- 0
        } else {
          fta = f(i/n-q)
          f_a[i] <- CDF(uniroot(fta, c(0,1), tol = tol)$root - q)
        }
        
        if((i-1)/n+q>1){
          f_b[i] <- 0
        } else {
          ftb = f((i-1)/n+q)
          f_b[i] <- CDF(uniroot(ftb, c(0,1), tol = tol)$root + q)
        }
      } 
    } else {
      for (i in 1:n) {
        if(i/n-q<0){
          f_a[i] <- 0
        } else {
          fta = f(i/n-q)
          temp = 0
          h = fta(temp)/pdf(temp)
          while (abs(h)>tol) {
            h = fta(temp)/pdf(temp)
            temp = temp - h
          }
          f_a[i] <- CDF(temp - q)
        }
        
        if((i-1)/n+q>1){
          f_b[i] <- 0
        } else {
          ftb = f((i-1)/n+q)
          temp = 0
          h = ftb(temp)/pdf(temp)
          while (abs(h)>tol) {
            h = ftb(temp)/pdf(temp)
            temp = temp - h
          }
          f_b[i] <- CDF(temp + q)
        }
      }
    }
    
    
  } else {
    for (i in 1:n) {
      f_b[i] <- CDF(CDFinverse((i-1)/n+q)+q)
      f_a[i] <- CDF(CDFinverse(i/n-q)-q)
    }
  }
  
  f_a[f_a<0] = 0
  f_b[f_b>1] = 1
  df <- data.frame(rbind(f_b, f_a))
  write.table(df, "Boundary_Crossing_Time.txt", 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)
  file.remove("Boundary_Crossing_Time.txt")
  
  names(q) <- "H"
  result = list(p.value = PVAL, method = "One-sample Hausdorff test",
                statistic = q, alternative = "two-sided",
                data.name = paste(deparse(substitute(x))))
  class(result) = "htest"
  return(result)
}

H_stat_1s_1d = function(x0, CDF, pdf = NULL, tol = 1e-10, max.init = 1000){
  
  f <- function(lambda){
    function(x){
      x + CDF(x) - lambda
    }
  }
  x0 = sort(x0)
  edfx = ecdf(x0)
  y0 = edfx(x0)
  yf = CDF(x0)
  n = length(x0)
  x0 = x0[ceiling(1:(2*n)/2)]
  y0 = c(0,y0[ceiling(1:(2*n-1)/2)])
  yf = yf[ceiling(1:(2*n)/2)]
  Index = rep(c(T,F),n)
  Index = (Index == (yf>y0))
  x0 = x0[Index]
  y0 = y0[Index]
  parameter = x0 + y0
  
  xs = rep(0,length(parameter))
  if(is.null(pdf)){
    for(i in 1:length(parameter)){
      ft = f(parameter[i])
      xs[i] = uniroot(ft, c(x0[i]-1,x0[i]+1), tol = tol)$root
    }
  } else {
    for(i in 1:length(parameter)){
      ft = f(parameter[i])
      h = ft(x0[i])/(1+pdf(x0[i]))
      xs[i] = x0[i]
      num.init = 0
      while (abs(h)>tol) {
        h = ft(xs[i])/(1+pdf(xs[i]))
        xs[i] = xs[i] - h
        num.init = num.init + 1
        if(max.init < num.init) stop("Maximum number of iteration reached.")
      }
    }
  }
  H_vec = xs - x0
  return(max(abs(H_vec)))
}



Hausdorff_test_1samp_simu = function(x, CDF, CDFinverse = NULL, nboots = 2000, tol = 1e-6){
  test.stat = H_stat_1s_1d(x, CDF, tol)
  n = length(x)
  bigger = 0L
  Inv = is.null(CDFinverse)
  if(!Inv){
    for(i in 1:nboots){
      x0 = runif(n)
      x0 = CDFinverse(x0)
      bigger = bigger + (H_stat_1s_1d(x0, CDF, tol)>test.stat)
    }
  } else {
    CDFinverse = function(x1) optim(0,(function (t) abs(CDF(t) - x1)), method="L-BFGS-B")$par
    for(i in 1:nboots){
      x0 = runif(n)
      x0 = unlist(lapply(x0,CDFinverse))
      bigger = bigger + (H_stat_1s_1d(x0, CDF, tol)>test.stat)
    }
  }
  
  names(test.stat) <- "H"
  result = list(p.value = bigger/nboots, method = "One-sample Hausdorff test",
                statistic = test.stat, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1))))
  class(result) = "htest"
  return(result)
}
