# Millstein version of Meff, effective number of independent variables
# rows of mydat are individuals, columns are variables
meff.jm = function(mydat, B = 1, seed){
  
  cmat = cor(mydat, use = "pairwise.complete.obs") 
  
  # remove missing
  arank1 = dim(cmat)[1]
  for(i in 1:arank1){
    if(dim(cmat)[2] >= i) { 
      navec = !is.na(cmat[,i])
      if(sum(as.integer(navec)) > 0) cmat = cmat[navec, navec]
    }
  }
  
  # compute eigenvalues
  rst = eigen(cmat, only.values = TRUE)
  v1 = rst$values
  
  set.seed(seed)
  out <- numeric()
  for(b in 1:B) {
    # randomize vector
    pfun = function(vec){
      vec = vec[order(runif(length(vec)))]
    }
    
    
    # Randomly permute each column wrt other columns enforcing independence
    mydat1 = apply(mydat,2,pfun)
    cmat = cor(mydat1, use = "pairwise.complete.obs") # rows of mydat are individuals
    
    # remove missing
    arank = dim(cmat)[1]
    for(i in 1:arank){
      if(dim(cmat)[2] >= i) { 
        navec = !is.na(cmat[,i])
        if(sum(as.integer(navec)) > 0) cmat = cmat[navec, navec]
      }
    }
    
    # compute eigen values for independent data
    rst = eigen(cmat, only.values = TRUE)
    v2 = rst$values
    
    # remove missing
    v1 = v1[!is.na(v1)]
    v2 = v2[!is.na(v2)]
    l = min(length(v1),length(v2))
    v1 = v1[1:l]
    v2 = v2[1:l]
    
    # variance inflation of top eigenvalues, sum1 vs sum2
    flag = 0; i = 1; sum1=0; sum2=0;
    while(flag == 0){
      if(v1[i] > v2[i]){
        sum1 = sum1 + v1[i]
        sum2 = sum2 + v2[i]
        if(i == l) flag = 1
      } else flag = 1
      i = i + 1
    }
    
    # effective no of independent vars is sum of eigenvalues for truly independent vars (n) 
    # minus the variance inflation of top eigenvalues, that is, (permuted - observed)
    # sum of eigen values of correlation matrix = trace = sum of diag = n
    # eigenvalue of variables X can be thought of as normalized % variance explained by X
    # dependencies between vars will increase top eigenvalues according to dependent %variance
    # therefore raw n = sum(v2) and dependent variance = (sum1 - sum2) or the inflation in eigenvalues due to dependence
    # thus, the effective number of independent vars = total variance - dependentent variance = sum(v2) - (sum1 - sum2)
    
    if(v1[1] > v2[1]){
      n = sum(v2) - (sum1 - sum2) 
    } else {
      n = arank1
    }
    
    out[b] <- n
  }
  
  return(mean(out))
}




