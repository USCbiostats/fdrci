#' Estimate the Effective Number of Tests
#' 
#' Estimate the effective number of tests using a permutation-based approach. 
#' 
#' @param mydat A design matrix with \eqn{n} observations (rows) and \eqn{p} covariates (columns).
#' @param B Integer. Number of permutations to perform. (Default is 1)
#' @param seed Integer. Setting the seed for reproducibility.
#' @return Numeric. Returns the estimated effective number of tests averaged over \code{B} permutations.
#' @author Joshua Millstein, Eric S. Kawaguchi
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#' @importFrom stats cor
#' @examples
#' 
#' # Independent
#' ss=100
#' nvar=100
#' X = as.data.frame(matrix(rnorm(ss * nvar), nrow = ss, ncol = nvar))
#' meff.jm(X, B = 5, seed = 1234)
#' 
#' # High correlation
#' S = matrix(0.9, ss, nvar)
#' diag(S) = 1
#' X = as.matrix(X) %*% chol(S)
#' meff.jm(X, B = 5, seed = 1234)
#' @export
#' 
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




