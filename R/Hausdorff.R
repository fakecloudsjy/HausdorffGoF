Hausdorff = function(a,b){
  edfx1 = ecdf(a)
  edfx2 = ecdf(b)
  a = unique(a)
  b = unique(b)
  xjoint = sort(c(a,b))
 # separate x1, x2 by observing the max and min of a and b
  ay1 = edfx1(xjoint)
  ay2 = edfx2(xjoint)
  tempy1 = 0
  tempy2 = 0
  x1 = x2 = y1 = y2 = NULL
  for(i in 1:length(c(a,b)) ){
    if(max(ay1[i],ay2[i]) > tempy1){
      tempy1 = max(ay1[i],ay2[i])
      x1 = c(x1, xjoint[i])
      y1 = c(y1,tempy1)
    } 
    if(min(ay1[i],ay2[i]) > tempy2){
      x2 = c(x2,xjoint[i])
      y2 = c(y2,tempy2)
      tempy2 = min(ay1[i],ay2[i])
    }
  }
  
  m1 = length(x1)
  m2 = length(x2)
  x1 = c(x1[ceiling((1:(2*m1))/2)], xjoint[length(xjoint)])
  y1 = c(0, y1[ceiling((1:(2*m1))/2)])
  x2 = c(xjoint[1],x2[ceiling((1:(2*m2))/2)])
  y2 = c(y2[ceiling((1:(2*m2))/2)],1)
  m1 = 2 * m1 + 1
  m2 = 2 * m2 + 1
  ##exchange vertices if former is larger than the other
  if(m1>m2){
    temp = x1
    x1 = - x2[m2:1]
    x2 = - temp[m1:1]
    temp = y1
    y1 = 1 - y2[m2:1]
    y2 = 1 - temp[m1:1]
    temp = m1
    m1 = m2
    m2 = temp
  }
  
  #give the parameter lambda 
  par1 = x1 + y1
  par2 = x2 + y2 # containing lambda_0
  par1 = par1[2*1:(m1%/%2)]
  ## par2 = 1:16
  ## par1 = c(2+1:3/3,3,4,4,9.5,10+1:6/6.1) 
  
  #calculate batch
  m = match(par1, sort(c(par1,par2)))
  temp = as.vector(table(m))
  temp2 = which(temp>1)
  Logic1 = length(temp2)>0
  if(Logic1){
    temp3 = cumsum(temp)
    for (i in temp2) {
      if(i == 1){
        m[1:temp2[1]] = 1:temp2[1]
      } else {
        m[(temp3[i-1]+1):temp3[i]] = 0:(temp3[i]-temp3[i-1]-1) + m[temp3[i]]
      }
      
    }
  }
  temp = which(diff(m)>1)
  n = c(temp,length(m)) #n
  ibatch = m[n] - n #i
  m = c(0,temp)+1 #m
  
  H_vec = rep(0,length(ibatch))
  #calculate Hausdorff distance in each batch
  for (i in 1:length(ibatch)) H_vec[i] = ifelse(ibatch[i]%%2,
                                                y1[2*n[i]] - y2[ibatch[i]+1],
                                                x2[ibatch[i]] - x1[2*m[i]]
                                                )
  return(max(H_vec))
}

Hausdorff_Rank = function(a,b){
  N = length(a) + length(b)
  xjoint = sort(c(a,b))
  a = match(a,xjoint)/N
  b = match(b,xjoint)/N
  xjoint = 1:N/N
  edfx1 = ecdf(a)
  edfx2 = ecdf(b)
  a = unique(a)
  b = unique(b)
  
  # separate x1, x2 by observing the max and min of a and b
  ay1 = edfx1(xjoint)
  ay2 = edfx2(xjoint)
  tempy1 = 0
  tempy2 = 0
  x1 = x2 = y1 = y2 = NULL
  for(i in 1:length(c(a,b)) ){
    if(max(ay1[i],ay2[i]) > tempy1){
      tempy1 = max(ay1[i],ay2[i])
      x1 = c(x1, xjoint[i])
      y1 = c(y1,tempy1)
    } 
    if(min(ay1[i],ay2[i]) > tempy2){
      x2 = c(x2,xjoint[i])
      y2 = c(y2,tempy2)
      tempy2 = min(ay1[i],ay2[i])
    }
  }
  
  m1 = length(x1)
  m2 = length(x2)
  x1 = c(x1[ceiling((1:(2*m1))/2)], xjoint[length(xjoint)])
  y1 = c(0, y1[ceiling((1:(2*m1))/2)])
  x2 = c(xjoint[1],x2[ceiling((1:(2*m2))/2)])
  y2 = c(y2[ceiling((1:(2*m2))/2)],1)
  m1 = 2 * m1 + 1
  m2 = 2 * m2 + 1
  ##exchange vertices if former is larger than the other
  if(m1>m2){
    temp = x1
    x1 = - x2[m2:1]
    x2 = - temp[m1:1]
    temp = y1
    y1 = 1 - y2[m2:1]
    y2 = 1 - temp[m1:1]
    temp = m1
    m1 = m2
    m2 = temp
  }
  
  #give the parameter lambda 
  par1 = x1 + y1
  par2 = x2 + y2 # containing lambda_0
  par1 = par1[2*1:(m1%/%2)]
  ## par2 = 1:16
  ## par1 = c(2+1:3/3,3,4,4,9.5,10+1:6/6.1) 
  
  #calculate batch
  m = match(par1, sort(c(par1,par2)))
  temp = as.vector(table(m))
  temp2 = which(temp>1)
  Logic1 = length(temp2)>0
  if(Logic1){
    temp3 = cumsum(temp)
    for (i in temp2) {
      if(i == 1){
        m[1:temp2[1]] = 1:temp2[1]
      } else {
        m[(temp3[i-1]+1):temp3[i]] = 0:(temp3[i]-temp3[i-1]-1) + m[temp3[i]]
      }
      
    }
  }
  temp = which(diff(m)>1)
  n = c(temp,length(m)) #n
  ibatch = m[n] - n #i
  m = c(0,temp)+1 #m
  
  H_vec = rep(0,length(ibatch))
  #calculate Hausdorff distance in each batch
  for (i in 1:length(ibatch)) H_vec[i] = ifelse(ibatch[i]%%2,
                                                y1[2*n[i]] - y2[ibatch[i]+1],
                                                x2[ibatch[i]] - x1[2*m[i]]
  )
  return(max(H_vec))
}





KS = function(a,b){
  edfx1 = ecdf(a)
  edfx2 = ecdf(b)
  xjoint = sort(c(a,b))
  return(max(abs(edfx1(xjoint)-edfx2(xjoint))))
}

Hausdorff_geq = function(a,b,q){ #examine whether H(a,b)>q
  edfx1 = ecdf(a)
  edfx2 = ecdf(b)
  a = unique(a)
  b = unique(b)
  xjoint = sort(c(a,b))
  # separate x1, x2 by observing the max and min of a and b
  ay1 = edfx1(xjoint)
  ay2 = edfx2(xjoint)
  tempy1 = 0
  tempy2 = 0
  x1 = x2 = y1 = y2 = NULL
  for(i in 1:length(c(a,b)) ){
    if(max(ay1[i],ay2[i]) > tempy1){
      tempy1 = max(ay1[i],ay2[i])
      x1 = c(x1, xjoint[i])
      y1 = c(y1,tempy1)
    } 
    if(min(ay1[i],ay2[i]) > tempy2){
      x2 = c(x2,xjoint[i])
      y2 = c(y2,tempy2)
      tempy2 = min(ay1[i],ay2[i])
    }
  }
  
  m1 = length(x1)
  m2 = length(x2)
  x1 = c(x1[ceiling((1:(2*m1))/2)], xjoint[length(xjoint)])
  y1 = c(0, y1[ceiling((1:(2*m1))/2)])
  x2 = c(xjoint[1],x2[ceiling((1:(2*m2))/2)])
  y2 = c(y2[ceiling((1:(2*m2))/2)],1)
  m1 = 2 * m1 + 1
  m2 = 2 * m2 + 1
  ##exchange vertices if former is larger than the other
  if(m1>m2){
    temp = x1
    x1 = - x2[m2:1]
    x2 = - temp[m1:1]
    temp = y1
    y1 = 1 - y2[m2:1]
    y2 = 1 - temp[m1:1]
    temp = m1
    m1 = m2
    m2 = temp
  }
  
  #give the parameter lambda 
  par1 = x1 + y1
  par2 = x2 + y2 # containing lambda_0
  par1 = par1[2*1:(m1%/%2)]
  ## par2 = 1:16
  ## par1 = c(2+1:3/3,3,4,4,9.5,10+1:6/6.1) 
  
  #calculate batch
  m = match(par1, sort(c(par1,par2)))
  temp = as.vector(table(m))
  temp2 = which(temp>1)
  Logic1 = length(temp2)>0
  if(Logic1){
    temp3 = cumsum(temp)
    for (i in temp2) {
      if(i == 1){
        m[1:temp2[1]] = 1:temp2[1]
      } else {
        m[(temp3[i-1]+1):temp3[i]] = 0:(temp3[i]-temp3[i-1]-1) + m[temp3[i]]
      }
      
    }
  }
  temp = which(diff(m)>1)
  n = c(temp,length(m)) #n
  ibatch = m[n] - n #i
  m = c(0,temp)+1 #m
  
  #calculate Hausdorff distance in each batch
  for (i in 1:length(ibatch)){
    temp = ifelse(ibatch[i]%%2,
                  y1[2*n[i]] - y2[ibatch[i]+1],
                  x2[ibatch[i]] - x1[2*m[i]]
    )
    if(temp > q) return(T)
  }
  return(F)
}



