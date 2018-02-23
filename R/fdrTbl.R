#' FDR Estimate and Confidence Interval Sequence Table
#' 
#' Computes FDR estimates and confidence intervals for a sequence of potential
#' significance thresholds.
#' 
#' fdrTbl calls fdr_od. Output from fdrTbl() can be used for FDRplot() input.
#' 
#' @param obs_vec observed vector of p-values.
#' @param perm_list list of dataframes that include a column of permutation
#' p-values (or statistics) in each. The length of the list permp = number of
#' permutations.
#' @param pname name of column in each list component dataframe that includes
#' p-values (or statistics).
#' @param ntests total number of observed tests, which is usually the same as
#' the length of obs_vec and the number of rows in each perm_list dataframe.
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
#' @param c1 overdispersion parameter. If this parameter is not specified
#' (default initial value is NA), then the parameter is estimated from the
#' data. If all tests are known to be independent, then this parameter should
#' be set to 1.
#' @return A dataframe is returned where rows correspond to p-value thresholds
#' in the sequence from lowerbound to upperbound and columns are:
#' c("threshold","fdr","ll","ul","pi0","odp","S","Sp") \item{threshold
#' }{p-value threshold chosen to define positive tests} \item{fdr }{estimated
#' FDR at the chosen p-value threshold} \item{ll }{estimated lower 95\%
#' confidence bound for the FDR estimate} \item{ul }{estimated upper 95\%
#' confidence bound for the FDR estimate} \item{pi0 }{estimated percent of true
#' null hypotheses} \item{odp }{estimated over-dispersion parameter} \item{S
#' }{observed number of positive tests} \item{Sp }{total number of positive
#' tests summed across all permuted result sets}
#' @author Joshua Millstein
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#' @keywords htest nonparametric
#' @examples
#' 
#' 
#' nrow_=100
#' ncol_=100
#' X = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
#' Y = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
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
#' for(p_ in 1:nperm){
#' 	X1 = X[order(runif(ncol_)),]
#' 	perml[[p_]] = myanalysis(X1,Y)
#' }
#' 
#' ## FDR results table
#' fdrTbl(obs$pvalue,perml,"pvalue",ncol_,1,2)
#' 
#' @export
fdrTbl <-
function(obs_vec,perm_list,pname,ntests,lowerbound,upperbound,incr=.1,cl=.95,c1=NA){
	# obs_vec is a vector of observed p-values
	# lowerbound and upperbound define -log10(p-value) range over which fdr is computed for a sequence of thresholds
	# If obs_vec and perm_list have high p-values filtered out, then lowerbound should be >= filtering threshold on -log10(p-value) scale 
	# perm_list is a list of dataframes, each of which include a column with permutation p-values
	# pname contains a string that is the name of the permutation p-value column
	# ntests is the number of tests conducted (if no filtering is done, ntests == length(obs_vec)

	plotdat = as.data.frame(matrix(NA,nrow=0,ncol=8))
	names(plotdat) = c("threshold","fdr","ll","ul","pi0","odp","S","Sp")
	thres_ = seq(lowerbound,upperbound,incr)
	for(i in 1:length(thres_)){
	   plotdat[i,"threshold"] = thres_[i]
	   thr = 10^-(thres_[i])
	   tmp = fdr_od(obs_vec,perm_list,pname,ntests,thr,cl=.95,c1)
	   if(length(tmp) == length(plotdat) - 1) plotdat[i,-1] = tmp
	}
	
	return(plotdat)
} # End fdrTbl
