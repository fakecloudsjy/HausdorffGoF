Hausdorff45 = function(a,b){
  m = length(a)
  n = length(b)
  ay1 = as.numeric(table(a))/m
  by1 = as.numeric(table(b))/n
  a = sort(unique(a))
  b = sort(unique(b))
  if(length(a)>length(b)){
    temp = a
    a = b
    b = temp
    temp = ay1
    ay1 = by1
    by1 = temp
  }
  min_val = min(c(a,b))
  a = a - min_val
  b = b - min_val
  max_val = max(c(a,b))
  a_len = length(a)
  b_len = length(b)
  aa = rep(0,2*a_len)
  bb = rep(0,2*b_len)
  
  aa[2*1:a_len] = ay1
  aa[2*1:a_len-1] = -c(a[1],diff(a))
  bb[2*1:b_len] = by1
  bb[2*1:b_len-1] = -c(b[1],diff(b))
  m = 2 * a_len
  if(a[a_len]==max_val){
    aa = aa[1:(2*a_len-1)]
    m = m - 1
    if(bb[2*b_len-1]!=-max_val) bb = c(bb,max_val)
  }

  
  ax = cumsum(abs(aa))
  ay = cumsum(aa)
  bx = cumsum(abs(bb))
  by = cumsum(bb)
  b_len = length(bb)
  

  H_vec = NULL
  j = 1
  begining = min(2, m)
  if(ay[1] != 0)  begining = 1
  while (T) {
    if(j > b_len) break
    if(bx[j]>ax[begining]) break
    j = j + 1
  }
  
  for(i in begining:m){

    H_temp = ay[i] - by[max(j-1,1)] + ifelse(j%%2,1,-1)*(ax[i] - bx[max(j-1,1)])
    if((H_temp<0) == (i%%2)) H_vec = c(H_vec, H_temp)
      
      if(i < m){
        while (T) {
          if(j > b_len) break
          if(bx[j]>ax[i+1]) break
          j = j + 1
        }
      }
  }

  out = c(max(abs(H_vec))/2
          )

  return(out)
}


