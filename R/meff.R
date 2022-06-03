#' Estimate the Effective Number of Tests
#' 
#' Estimate the effective number of tests using a permutation-based approach. 
#' 
#' The effective no of independent vars is the sum of the scaled eigenvalues for truly independent vars (n) 
#' minus the variance inflation of the top eigenvalues, that is, (observed - permuted)
#' The sum of the eigenvalues of the correlation matrix = trace = sum of diag = n.
#' The eigenvalues of variables X can be thought of as normalized % variance explained by X.
#' Dependencies between variables will increase the magnitude of top eigenvalues,
#' therefore, the effective number of independent variables is equal to the proportion of the
#' sum of eigenvalues attributed to independence minus the proportion attributed to dependencies.
#' 
#' @param mydat A design matrix with \eqn{n} observations (rows) and \eqn{p} covariates (columns).
#' @param B Integer. Number of permutations to perform. (Default is 1)
#' @return Numeric. Returns the estimated effective number of tests averaged over \code{B} permutations.
#' @author Joshua Millstein, Eric S. Kawaguchi
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#' @importFrom stats cor runif
#' @examples
#' 
#' # Independent
#' ss=300
#' nvar=100
#' X = as.data.frame(matrix(rnorm(ss * nvar), nrow = ss, ncol = nvar))
#' meff.jm(X, B = 5)
#' 
#' # High correlation
#' S = matrix(0.9, nvar, nvar)
#' diag(S) = 1
#' X = as.matrix(X) %*% chol(S)
#' meff.jm(X, B = 5)
#' @export
#' 
meff.jm = function(mydat, B = 1){
  
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
    
    if(v1[1] > v2[1]){
      n = sum(v2) - (sum1 - sum2) 
    } else {
      n = arank1
    }
    
    out[b] <- n
  }
  
  return(mean(out))
}





