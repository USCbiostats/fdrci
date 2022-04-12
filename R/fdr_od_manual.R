#' Tail-area FDR and Confidence Interval (User customizable)
#' 
#' This function can be used to estimate FDR, corresponding confidence
#' interval for a given estimate of the number of tests and number of rejected huptheses.
#' 
#' If a very large number of tests are conducted, it may be useful to filter
#' results, that is, save only results of those tests that meet some relaxed
#' nominal significance threshold. This alleviates the need to record results
#' for tests that are clearly non-significant. Results from fdr_od() are valid
#' as long as thres < the relaxed nomimal significance threshold for both
#' observed and permuted results. It is not necessary for the input to fdr_od()
#' to be p-values, however, fdr_od() is designed for statistics in which
#' smaller values are more extreme than larger values as is the case for
#' p-values. Therefore, if raw statistics are used, then a transformation may
#' be necessary to insure that smaller values are more likely associated with
#' false null hypotheses than larger values. In certain situations, for
#' instance when a large proportion of tests meet the significance threshold,
#' \code{pi0} is estimated to be very small, and thus has a large influence on the FDR
#' estimate. To limit this influence, \code{pi0} is constrained to be .5 or greater,
#' resulting in a more conservative estimate under these conditions.
#' 
#' @param obsp observed vector of p-values.
#' @param thres significance threshold.
#' @param cl confidence level (default is .95).
#' @param m Numeric. User-provided estimate of the number of the tests. Default is \code{NULL}, which implies the standard calculation (\code{m = length(obsp)}).
#' @param s Numeric. User-provided estimate of the number of rejected hypotheses. Default is \code{NULL}, which implies the standard calculation (\code{s = sum(obsp <= thres)}).
#' @return 
#' A list which includes: \item{FDR }{FDR point estimate} \item{ll
#' }{lower confidence limit} \item{ul }{upper confidence limit} \item{M}{estimated total number of tests}
#' \item{S }{estimated number of positive tests (rejected hypotheses)}
#' 
#' @author Joshua Millstein, Eric S. Kawaguchi
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#' @importFrom stats qnorm var
#' @keywords htest nonparametric
#' @examples
#' 
#' ss=100
#' nvar=100
#' X = as.data.frame(matrix(rnorm(ss*nvar),nrow=ss,ncol=nvar))
#' e = as.data.frame(matrix(rnorm(ss*nvar),nrow=ss,ncol=nvar))
#' Y = .1*X + e
#' nperm = 10
#' 
#' myanalysis = function(X,Y){
#' 	ntests = ncol(X)
#' 	rslts = as.data.frame(matrix(NA,nrow=ntests,ncol=2))
#' 	names(rslts) = c("ID","pvalue")
#' 	rslts[,"ID"] = 1:ntests
#' 	for(i in 1:ntests){
#' 		fit = cor.test(X[,i],Y[,i],na.action="na.exclude",
#' 			alternative="two.sided",method="pearson")
#' 		rslts[i,"pvalue"] = fit$p.value
#' 	}
#' 	return(rslts)
#' } # End myanalysis
#' 
#' # Generate observed results
#' obs = myanalysis(X,Y)
#' 
#' fdr_od_manual(obs$pvalue, thres = .05, cl = 0.95, m = NULL, s = NULL)
#' 
#' # User-inputted m and s
#' fdr_od_manual(obs$pvalue, thres = .05, cl = 0.95, m = 100, s = 10)
#' 
#' @export
fdr_od_manual <-
  function(obsp, 
           thres = 0.05,
           cl = .95,
           m = NULL,
           s = NULL){
    
    # Checks
    if(cl < 0 | cl > 1) stop("cl must be between 0 and 1")
    if(thres < 0 | thres > 1) stop("thres must be between 0 and 1")
    
    z_ = qnorm(1 - (1 - cl)/2) # two-tailed test
    
    s0 = sum(obsp <= thres)
    m0 = length(obsp)
    
    if(is.null(m)) {
      m = m0
    }
    
    if(is.null(s)) {
      s = s0
    }
    
    if(s0 < s) stop("User-inputted value for s must be less than or equal to the observed number of rejected hypotheses")
    if(m0 < m) stop("User-inputted value for m must be less than or equal to the total number of tests")
    if(m < 1) stop("m must be NULL or greater than or equal to 1")
    if(s < 1) stop("s must be NULL or greater than or equal to 1")
    if(m < s) stop("m must be greater than or equal to s")
    
    # If thres is 1 then FDR is 1
    if (thres == 1) {
        fdr = 1
      } else {
        t1 = m * thres / s
        t2 = (1 - (s / m)) / (1 - thres)
        fdr = t1 * t2
    }
      
      if(fdr > 1) fdr = 1
      
      s2fdr = (1 / s) + (1 / (m - s))
      
      # Get upper and lower limit
      ul = exp(log(fdr) + z_ * sqrt(s2fdr))
      ll = exp(log(fdr) - z_ * sqrt(s2fdr))
      if(ul > 1) ul = 1
      
      rslt = c(fdr, ll, ul, m, s)
      names(rslt) = c("fdr", "ll", "ul", "M", "S")
      
      return(rslt)
  }