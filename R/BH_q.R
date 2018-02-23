#' Generate q-values based on the BH (Benjamini and Hochberg) approach
#' 
#' Benjamini and Hochberg (1995) developed the first approach for controlling FDR. Unlike the approach proposed by Millstein and Volfson (2013), it is based on the idea that a significance level (alpha) is specified a priori by the investigator. Storey and Tibshirani (2003) developed the q-value concept, where FDR is estimated at each observed p-value.
#' 
#' @param pvals vector of observed p-values.
#' @author Joshua Millstein
#' @references Benjamini, Yoav, and Yosef Hochberg. "Controlling the false discovery rate: a practical and powerful approach to multiple testing." Journal of the royal statistical society. Series B (Methodological) (1995): 289-300.
#' 
#' Storey, John D., and Robert Tibshirani. "Statistical significance for genomewide studies." Proceedings of the National Academy of Sciences 100.16 (2003): 9440-9445.
#'
#' 
#' @examples
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
#' # Generate observed results
#' obs = myanalysis(X,Y)
#' q.values.BH = BH_q( obs[, "pvalue"] )
#' 
#' @export
BH_q <-
function(pvals){            
	m = length(pvals)
	po = pvals[order(pvals)]
	i = 1:m
	new.order = i[order(pvals)]
	qval = po * m /  i
	for(i in 1:m){
		qval[ 1:i ] = ifelse( qval[ 1:i ] > qval[ i ], qval[ i ], qval[ 1:i ] )
		if( qval[ i ] > .9 ) break
	}
	ind = order( new.order )  # reorder back to original
	qval = qval[ ind ]
	return(qval)
} # End BH_q function

