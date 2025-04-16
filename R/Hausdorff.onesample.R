
Hausdorff1sample = function(x0, CDF, pdf = NULL, tol = 1e-10){
  
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
      while (abs(h)>tol) {
        h = ft(xs[i])/(1+pdf(xs[i]))
        xs[i] = xs[i] - h
      }
    }
  }
  H_vec = xs - x0
  return(max(abs(H_vec)))
}

Hausdorff_test_1sample = function(x, CDF, CDFinverse = NULL, pdf = NULL, tol = 1e-10){
  q = Hausdorff1sample(x, CDF, pdf, tol)
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
                data.name = paste(deparse(substitute(x1))))
  class(result) = "htest"
  return(result)
}

Hausdorff_test_1samp_simu = function(x, CDF, CDFinverse = NULL, nboots = 2000, tol = 1e-6){
  test.stat = Hausdorff1sample(x, CDF, tol)
  n = length(x)
  bigger = 0L
  Inv = is.null(CDFinverse)
  if(!Inv){
    for(i in 1:nboots){
      x0 = runif(n)
      x0 = CDFinverse(x0)
      bigger = bigger + (Hausdorff1sample(x0, CDF, tol)>test.stat)
    }
  } else {
    CDFinverse = function(x1) optim(0,(function (t) abs(CDF(t) - x1)), method="L-BFGS-B")$par
    for(i in 1:nboots){
      x0 = runif(n)
      x0 = unlist(lapply(x0,CDFinverse))
      bigger = bigger + (Hausdorff1sample(x0, CDF, tol)>test.stat)
    }
  }
  
  names(test.stat) <- "H"
  result = list(p.value = bigger/nboots, method = "One-sample Hausdorff test",
                statistic = test.stat, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1))))
  class(result) = "htest"
  return(result)
}

Hausdorff1sample_R = function(x0, CDF, pdf = NULL, CDFinverse = NULL, tol = 1e-10){
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
  x0_y0inv = CDFinverse(y0)
  parameter = x0 + y0
  H_vec = rep(0,length(parameter))
  for(i in 1:length(parameter)){
    if(x0[i]==x0_y0inv[i]){
      H_vec[i] = 0
    } else {
      ft = function(x) y0[i] - CDF(x) - pdf(x)*(x - x0[i])
      if(y0[i] == 1){
        H_vec[i] = -optim(x0[i],(function (h) -(y0[i] - CDF(h))*(h - x0[i])), method="L-BFGS-B", lower = 0)$value
      } else {
        h = uniroot(ft, sort(c(x0[i],min(x0_y0inv[i],1e10))), tol = tol)$root
        H_vec[i] = (y0[i] - CDF(h))*(h - x0[i])
      }
    }
  }
  return(max(abs(H_vec)))
}

Hausdorff_Rec_1sample = function(x0, tol = 1e-10, pdf = NULL, CDFinverse = NULL, null, ...){
  if(is.character(null)){
    if(substr(null,1,1)=='p') y = substr(null,2,nchar(null))
    CDF <- get(paste("p",y,sep = ""), mode = "function", envir = parent.frame())
    pdf <- get(paste("d",y,sep = ""), mode = "function", envir = parent.frame())
    CDFinverse <- get(paste("q",y,sep = ""), mode = "function", envir = parent.frame())
  } else {
    CDF <- null
  }
  
  q = Hausdorff1sample_R(x0, CDF = CDF,pdf = pdf, CDFinverse = CDFinverse, tol = tol)
  n = length(x0)
  f_b = rep(0,n)
  f_a = rep(0,n)
  
  CDFinv = function(x) return(ifelse(x>1,Inf,ifelse(x<0,-Inf,CDFinverse(x))))
  for (i in 1:n) {
    fun1 = function(h){
      return(-CDF(CDFinv(i/n-sqrt(h*q)) - sqrt(q/h)))
    }
    fun2 = function(h){
      return(CDF(CDFinv((i-1)/n+sqrt(h*q))+sqrt(q/h)))
    }
    f_a[i] = -optim(x0[i], fun1, method="L-BFGS-B", lower = tol)$value
    f_b[i] = optim(x0[i], fun2, method="L-BFGS-B", lower = tol,upper = (1-(i-1)/n)^2/q)$value
  }
  
  
  f_a[f_a<0] = 0
  f_b[f_b>1] = 1
  df <- data.frame(rbind(f_b, f_a))
  write.table(df, "Boundary_Crossing_Time.txt", 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)
  file.remove("Boundary_Crossing_Time.txt")
  
  names(q) <- "H_Rec"
  result = list(p.value = PVAL, method = "One-sample Hausdorff-Rectangle test",
                statistic = q, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1))))
  class(result) = "htest"
  return(result)
}

Hausdorff_Rec_1sample2 = function(x0, tol = 1e-10, pdf = NULL, CDFinverse = NULL, null, ...){
  if(is.character(null)){
    if(substr(null,1,1)=='p') y = substr(null,2,nchar(null))
    CDF <- get(paste("p",y,sep = ""), mode = "function", envir = parent.frame())
    pdf <- get(paste("d",y,sep = ""), mode = "function", envir = parent.frame())
    CDFinverse <- get(paste("q",y,sep = ""), mode = "function", envir = parent.frame())
  } else {
    CDF <- null
  }
  
  q = Hausdorff1sample_R(x0, CDF = CDF,pdf = pdf, CDFinverse = CDFinverse, tol = tol)
  n = length(x0)
  f_b = rep(0,n)
  f_a = rep(0,n)
  
  CDFinv = function(x) return(ifelse(x>1,Inf,ifelse(x<0,-Inf,CDFinverse(x))))
  for (i in 1:n) {
    fun1 = function(h){
      return(pdf(CDFinv(i/n-sqrt(h*q)) - sqrt(q/h))-h)
    }
    fun2 = function(h){
      return(pdf(CDFinv((i-1)/n+sqrt(h*q))+sqrt(q/h))-h)
    }
    opt.h1 = uniroot(fun1,lower = (i/n)^2/q, tol = tol)
    opt.h2 = uniroot(fun2,upper = (1-(i-1)/n)^2/q, tol = tol)
    f_a[i] = CDF(CDFinv(i/n-sqrt(opt.h1*q)) - sqrt(q/opt.h1))
    f_b[i] = CDF(CDFinv((i-1)/n+sqrt(opt.h2*q))+sqrt(q/opt.h2))
  }
  
  
  f_a[f_a<0] = 0
  f_b[f_b>1] = 1
  df <- data.frame(rbind(f_b, f_a))
  write.table(df, "Boundary_Crossing_Time.txt", 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)
  file.remove("Boundary_Crossing_Time.txt")
  
  names(q) <- "H_Rec"
  result = list(p.value = PVAL, method = "One-sample Hausdorff-Rectangle test",
                statistic = q, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1))),
                non.crossing = df
                )
  class(result) = "htest"
  return(result)
}
