#' Permutation-Based FDR and Confidence Interval
#' 
#' This function can be used to estimate FDR, corresponding confidence
#' interval, and pi0, the proportion of true null hypotheses, given a selected
#' significance threshold, and results from permuted data.
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
#' pi0 is estimated to be very small, and thus has a large influence on the FDR
#' estimate. To limit this influence, pi0 is constrained to be .5 or greater,
#' resulting in a more conservative estimate under these conditions.
#' 
#' @param obsp observed vector of p-values.
#' @param permp list of dataframes that include a column of permutation
#' p-values (or statistics) in each. The length of the list permp = number of
#' permutations.
#' @param pnm name of column in each list component dataframe that includes
#' p-values (or statistics).
#' @param ntests total number of observed tests, which is usually the same as
#' the length of obsp and the number of rows in each permp dataframe. However,
#' this may not be the case if results were filtered by a p-value threshold or
#' statistic threshold. If filtering was conducted then thres must be smaller
#' (more extreme) than the filtering criterion.
#' @param thres significance threshold.
#' @param cl confidence level (default is .95).
#' @param c1 overdispersion parameter. If this parameter is not specified
#' (default initial value is NA), then the parameter is estimated from the
#' data. If all tests are known to be independent, then this parameter should
#' be set to 1.
#' @return A list which includes: \item{FDR }{FDR point estimate} \item{ll
#' }{lower confidence limit} \item{ul }{upper confidence limit} \item{pi0
#' }{proportion of true null hypotheses} \item{c1 }{overdispersion parameter}
#' \item{S }{observed number of positive tests} \item{Sp }{total number of
#' positive tests summed across all permuted result sets}
#' @author Joshua Millstein
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
#' ## Generate permuted results
#' perml = vector('list',nperm)
#' for(p_ in 1:nperm){
#' 	X1 = X[order(runif(nvar)),]
#' 	perml[[p_]] = myanalysis(X1,Y)
#' }
#' 
#' ## FDR results
#' fdr_od(obs$pvalue,perml,"pvalue",nvar,.05)
#' @export
fdr_od <-
function(obsp,permp,pnm,ntests,thres,cl=.95,c1=NA){
      z_ = qnorm(1 - (1 - cl)/2) # two-tailed test
	pcount = rep(NA,length(permp))
	for(p_ in 1:length(permp)){
	   permp[[p_]][,pnm] = ifelse(permp[[p_]][,pnm] <= thres,1,0)
	   pcount[p_] = sum(permp[[p_]][,pnm],na.rm=TRUE)
	}
	# over-dispersion parameter is observed variance of p in permuted data / expected
	p = mean(pcount,na.rm=TRUE)/ntests # estimate p
	e_vr = ntests*p*(1 - p)
	o_vr = var(pcount,na.rm=TRUE)
	if( is.na( c1 ) ) {
	   c1 = o_vr / e_vr
	   if( !is.na( c1 ) ) if( c1 < 1) c1 = 1
	} 

	nperm = length(permp)
	mo = ntests
	ro = sum(obsp <= thres)
	vp = sum(pcount)
	vp1 = vp
  	rslt = rep(NA,4)
  	if(ro > 0){
		if(vp == 0) vp = 1
		mean.vp = vp / nperm
        		fdr0 = mean.vp / ro
		num = mo - ro
		denom = mo - (vp/nperm)
		pi0 = num / denom
		
		# Rule of thumb: if num or denom is < 20, then pi0 is not stable so set it to the conservative value of 1
		pi0 = ifelse( num < 20 | denom < 20, 1, pi0 )
		if( pi0 > 1 ) pi0 = 1
        		if( pi0 < 0.5 ) pi0 = 0.5    # updated 10/14/16 to limit influence of pi0
        		fdr = fdr0 * pi0    # updated calculation of fdr to be robust to ro = mtests
		
		# variance of FDR
        		mp = nperm * mo
        		t1 = 1 / vp
        		denom = mp - vp
        		t2 = 1 / denom
        		t3 = 1 / ro
        		denom = ntests - ro
        		if( denom < 1 ) denom = 1  # updated 10/14/16 to avoid inf t4
        		t4 = 1 /  denom
        		s2fdr = (t1 + t2 + t3 + t4) * c1
        		ul = exp(log(fdr) + z_ * sqrt(s2fdr))
        		ll = exp(log(fdr) - z_ * sqrt(s2fdr))
	
		rslt = c(fdr,ll,ul,pi0)
		rslt = ifelse(rslt > 1, 1, rslt) # FDR > 1 does not make sense, thus set to 1 in this case
		rslt = c(rslt,c1,ro,vp1)
		names(rslt) = c("fdr", "ll", "ul", "pi0", "c1", "S", "Sp")
	} 
	return(rslt)
}

