#' FDR Estimate and Confidence Interval Sequence Table
#' 
#' Computes FDR estimates and confidence intervals for a sequence of potential
#' significance thresholds.
#' 
#' fdrTbl calls fdr_od for a series of discovery thresholds. 
#' Output from fdrTbl() can be used for FDRplot() input.
#' 
#' @param obs.vec observed vector of p-values.
#' @param perm.list list of dataframes that include a column of permutation
#' p-values (or statistics) in each. The length of the list permp = number of
#' permutations.
#' @param pname name of column in each list component dataframe that includes
#' p-values (or statistics).
#' @param ntests total number of observed tests, which is usually the same as
#' the length of obs.vec and the number of rows in each perm.list dataframe.
#' However, this may not be the case if results were filtered by a p-value
#' threshold or statistic threshold. If filtering was conducted then lowerbound
#' must be greater (more extreme) than the filtering criterion.
#' @param lowerbound lowerbound refers to the range of -log10(p-value) over
#' which fdr is computed for a sequence of thresholds
#' @param upperbound upperbound refers to the range of -log10(p-value) over
#' which fdr is computed for a sequence of thresholds
#' @param incr value by which to increment the sequence from lowerbound to
#' upperbound on a -log10(p-value) scale. Default is 0.1.
#' @param cl confidence level (default is .95).
#' @param c1 overdispersion parameter to account for dependencies among tests. If 
#' all tests are known to be independent, then this parameter should be set to 1.
#' @param meff (For parametric estimation, if \code{perm.list = NULL}.) Logical. To be passed into \code{fdr_od}. \code{TRUE} implies the calculation of the effective number of tests based on the JM estimator (Default is \code{TRUE})
#' @param seff (For parametric estimation, if \code{perm.list = NULL}.) Logical. To be passed into \code{fdr_od}. \code{TRUE} implies the calculation of the effective number of rejected hypotheses based on the JM estimator (Default is \code{TRUE})
#' @param mymat (For parametric estimation, if \code{perm.list = NULL}.) Matrix. To be passed into \code{fdr_od}. Design matrix used to calculate the p-values provided in \code{obsp}.
#' @param nperms (For parametric estimation, if \code{perm.list = NULL}.) Integer. To be passed into \code{fdr_od}. Number of permutations needed to estimate the effective number of (rejected) tests. (Must be non-zero, default is 5)
#' @param correct {"none", "BH"}, should confidence intervals be corrected for 
#' multiplicity using a modification of the Benjamini and Yekutieli (2005) approach 
#' for selecting and correcting intervals? (default is "none")
#' @details If correct = "BH", then confidence intervals will be corrected according to
#' the thresholds specified by lowerbound, upperbound, and incr. Thresholds will 
#' be selected if FDR is determined to be significantly different than 1. First
#' a Z-score test is conducted using the Millstein & Volfson standard error estimate.
#' Then BH FDR is computed according to the Benjamini and Yekutieli (2005) approach.
#' CIs for selected thresholds will be adjusted to account for multiple CI estimation.
#' For thresholds that are not selected, NA values are returned.
#' @return A dataframe is returned where rows correspond to p-value thresholds
#' in the sequence from lowerbound to upperbound and columns are:
#' 
#' If permutation:
#' c("threshold","fdr","ll","ul","pi0","odp","S","Sp") \item{threshold
#' }{p-value threshold chosen to define positive tests} \item{fdr }{estimated
#' FDR at the chosen p-value threshold} \item{ll }{estimated lower 95\%
#' confidence bound for the FDR estimate} \item{ul }{estimated upper 95\%
#' confidence bound for the FDR estimate} \item{pi0 }{estimated percent of true
#' null hypotheses} \item{odp }{estimated over-dispersion parameter} \item{S
#' }{observed number of positive tests} \item{Sp }{total number of positive
#' tests summed across all permuted result sets}
#' 
#' If parametric:
#' c("threshold","fdr","ll","ul","M","M.eff","S","S.eff") \item{threshold
#' }{p-value threshold chosen to define positive tests} \item{fdr }{estimated
#' FDR at the chosen p-value threshold} \item{ll }{estimated lower 95\%
#' confidence bound for the FDR estimate} \item{ul }{estimated upper 95\%
#' confidence bound for the FDR estimate} \item{M }{total number of tests}
#' \item{M.eff}{Effective number of tests via the JM estimator} \item{S
#' }{observed number of positive tests} \item{S.eff }{effective number of positive tests via the JM estimator}
#' 
#' @importFrom stats p.adjust
#' @author Joshua Millstein, Eric S. Kawaguchi
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#'
#' Benjamini, Yoav, and Daniel Yekutieli. "False discovery rate adjusted multiple 
#' confidence intervals for selected parameters." Journal of the American Statistical 
#' Association 100.469 (2005): 71-81.
#' 
#' @keywords htest nonparametric
#' @examples
#' 
#' 
#' n.row=100
#' n.col=100
#' X = as.data.frame(matrix(rnorm(n.row*n.col),nrow=n.row,ncol=n.col))
#' e = as.data.frame(matrix(rnorm(n.row*n.col),nrow=n.row,ncol=n.col))
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
#' ## Generate observed results
#' obs = myanalysis(X,Y)
#' 
#' ## Generate permuted results
#' perml = vector('list',nperm)
#' for(perm in 1:nperm){
#' 	X1 = X[order(runif(n.col)),]
#' 	perml[[perm]] = myanalysis(X1,Y)
#' }
#' 
#' ## FDR results table
#' fdrTbl(obs$pvalue,perml,"pvalue",n.col,1,2)
#' fdrTbl(obs$pvalue,perml,"pvalue",n.col,1,2,correct="BH")
#' 
#' ## FDR results table (parametric)
#' fdrTbl(obs$pvalue, NULL, "pvalue",n.col,1,2,meff = TRUE, seff = TRUE, mymat = X, nperms = 5)
#' 
#' @importFrom stats pnorm
#' @export
fdrTbl <-
function(obs.vec, perm.list = NULL, pname, ntests, lowerbound, upperbound, incr = .1, cl = .95, c1 = NA, correct = "none",
        # Options for parametric
        meff = TRUE, seff = TRUE, mymat, nperms = 5) {
    
	# obs.vec is a vector of observed p-values
	# lowerbound and upperbound define -log10(p-value) range over which fdr is computed for a sequence of thresholds
	# If obs.vec and perm.list have high p-values filtered out, then lowerbound should be >= filtering threshold on -log10(p-value) scale 
	# perm.list is a list of dataframes, each of which include a column with permutation p-values
	# pname contains a string that is the name of the permutation p-value column
	# ntests is the number of tests conducted (if no filtering is done, ntests == length(obs.vec)

  plotdat = as.data.frame(matrix(NA, nrow = 0, ncol = 8))
  thres = seq(lowerbound,upperbound,incr)
  
  if(is.null(perm.list)) {
    # Parametric 
    names(plotdat) = c("threshold", "fdr", "ll", "ul", "M", "M.eff", "S", "S.eff")
    
    # Checks
    if(nperms <= 0) stop("nperms must be a positive integer")
    if(!(meff %in% c(FALSE, TRUE))) stop("meff must be TRUE or FALSE")
    if(!(seff %in% c(FALSE, TRUE))) stop("seff must be TRUE or FALSE")
    if(!is.matrix(mymat) & !is.data.frame(mymat)) stop("mymat must be a matrix or data frame")
    
    for(i in 1:length(thres)) {
      plotdat[i,"threshold"] = thres[i]
      thr = 10^-(thres[i])
      tmp = fdr_od(obsp=obs.vec, permp = perm.list, pnm=pname, ntests=ntests, thres=thr, cl = .95, c1=c1, 
                   meff = meff, seff = seff, mymat = mymat, nperms = nperms)
      if(length(tmp) == length(plotdat) - 1) plotdat[i, -1] = tmp
    }
  } else {
    # Permutation
    names(plotdat) = c("threshold", "fdr", "ll", "ul", "pi0", "odp", "S", "Sp")
	
  	for(i in 1:length(thres)){
  	   plotdat[i,"threshold"] = thres[i]
  	   thr = 10^-(thres[i])
  	   tmp = fdr_od(obsp=obs.vec, permp = perm.list, pnm=pname, ntests=ntests, thres=thr, cl = .95, c1=c1)
  	   if(length(tmp) == length(plotdat) - 1) plotdat[i,-1] = tmp
  	}
  }
	
	if(is.element(correct, "BH")){
		aa = !is.na(plotdat$fdr)
		bb = !is.na(plotdat$ll)
		cc = !is.na(plotdat$ul)
		indx = which( aa & bb & cc )
		plotdat[!is.element(1:nrow(plotdat),indx), "ll"] = NA
    plotdat[!is.element(1:nrow(plotdat),indx), "ul"] = NA
    if (length(indx) > 1) {
      tmpt = plotdat[indx, c("fdr", "ll", "ul")]
      alpha = 1 - cl
      se.fdr = (log(tmpt$fdr) - log(tmpt$ll))/qnorm(1 - alpha/2)
      z = log(tmpt$fdr)/se.fdr
      p = pnorm(z, lower.tail = TRUE)
      qval = p.adjust(p, method="BH")
      sig = qval <= alpha
      R = sum(sig)
      alpha.a = R * alpha/length(p)        
      ll.a = exp(log(tmpt$fdr) - qnorm(1 - alpha.a/2) * se.fdr)
      ul.a = exp(log(tmpt$fdr) + qnorm(1 - alpha.a/2) * se.fdr)
      ll.a = ifelse(sig, ll.a, NA)
      ul.a = ifelse(sig, ul.a, NA)
      plotdat[indx, "ll"] = ifelse(ll.a <= 1, ll.a, 1)
      plotdat[indx, "ul"] = ifelse(ul.a <= 1, ul.a, 1)
		} # end if length indx
	} # end correct
	
	return(plotdat)
} # End fdrTbl


