
H_test_2s_1d = function(x1, x2, nboots = 2000, Exact = FALSE){
  test_stat = H_stat_2s_1d_tr(x1, x2) #Finds test stat
  comb = c(x1, x2)
  n = length(comb)
  na = min(length(x1),length(x2))
  if(choose(n,na)>2.147483646e9){
    Exact = F
    message("Sample sizes too large, switch to Monte Carlo")
  }
  if(!Exact){
    Method = "Monte Carlo"
    nboots = as.integer(nboots)					#Speeds up comparison below.
    reps = 0L
    bigger = 0L							  #Initializes Counter
    while (reps < nboots) {						#Loops over vector
      vec_labels = rep(F,n)
      vec_labels[sample.int(n,na,F)] = T #Samples indexes
      boot_t = H_stat_2s_1d_tr(comb[vec_labels],comb[!vec_labels]) #boot strap test stat
      bigger = bigger +(boot_t >= test_stat) #if new stat is bigger, increment
      reps = 1L+reps
    }
    out = bigger/nboots
    out[which(out==0)] = 1/(2*nboots)
    
  } else {
    Method = "Exact"
    bigger = 0L
    combn_exact = combn(n, na)
    for(i in 1:choose(n,na)){
      bigger = bigger + (H_stat_2s_1d_tr(comb[combn_exact[,i]], comb[-combn_exact[,i]]) >= test_stat)
    }
    out = bigger/choose(n,na)
    
  }
  names(test_stat) <- "H"
  result = list(p.value = out, method = paste("Two-sample Hausdorff Test (",Method, ")", sep = ""),
                statistic = test_stat, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1)), "and", deparse(substitute(x2))))
  class(result) = "htest"
  return(result)
}