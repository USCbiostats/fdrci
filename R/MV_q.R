#' MV q-values and confidence intervals
#' 
#' q-values with confidence intervals are generated, based in the Millstein and Volfson (MV) estimators. 
#'
#' @details 
#' Millstein and Volfson (2013) FDR is based on the idea that FDR is estimated at a level specified by the investigator. 
#' Storey and Tibshirani (2003) developed the q-value concept, where FDR is estimated at each observed p-value. 
#' However, Millstein and Volfson argued that in order to be informative, uncertainty in the estimate should be 
#' quantified, thus the development of confidence intervals for FDR. The MV FDR estimator is less conservative than the BH estimator.
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
#' @param cl confidence level (default is .95).
#' @param c1 overdispersion parameter. If this parameter is not specified
#' (default initial value is NA), then the parameter is estimated from the
#' data. If all tests are known to be independent, then this parameter should
#' be set to 1.
#' @return A dataframe which includes: \item{q}{q-value corresponding to the respective p-value} \item{q.ll}{q-value lower limit} \item{q.ul}{q-value upper limit} 
#' @author Joshua Millstein
#' @references Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. 
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#'
#' Storey, John D., and Robert Tibshirani. "Statistical significance for genomewide studies." 
#' Proceedings of the National Academy of Sciences 100.16 (2003): 9440-9445.
#' 
#' @examples
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
#' # Generate permuted results
#' perml = vector('list',nperm)
#' for(p_ in 1:nperm){
#' 	X1 = X[order(runif(nvar)),]
#' 	perml[[p_]] = myanalysis(X1,Y)
#' }
#' 
#' q.values.MV = MV_q(obs$pvalue,perml,"pvalue",nvar)
#' 
#' @export
MV_q <-
function(obsp,permp,pnm,ntests,cl=.95,c1=NA){            
	nms = c("p", "q", "q.ll", "q.ul")
	qvals = as.data.frame(matrix(NA, nrow=length(obsp), ncol=length(nms)))
	names(qvals) = nms
	qvals[, "p"] = obsp
	for(j in 1:length(obsp)){
		qvals[ j, c("q", "q.ll", "q.ul")] = fdr_od(obsp,permp,pnm,ntests,obsp[ j ])[c("fdr", "ll", "ul")]
	}
	# enforce monotonicity 
	op = order(obsp)
	for(tst in 1:nrow(qvals)){
		aa = qvals[ op[1:tst], "q" ] > qvals[ op[tst], "q" ]
		qvals[ op[1:tst], "q" ] = ifelse( aa, qvals[ op[tst], "q" ], qvals[ op[1:tst], "q" ] )
		qvals[ op[1:tst], "q.ll" ] = ifelse( aa, qvals[ op[tst], "q.ll" ], qvals[ op[1:tst], "q.ll" ] )
		qvals[ op[1:tst], "q.ul" ] = ifelse( aa, qvals[ op[tst], "q.ul" ], qvals[ op[1:tst], "q.ul" ] )
	}	
	return(qvals)
} # End MV_q function




