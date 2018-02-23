

#' Permutation-Based FDR Point and Confidence Interval Estimation
#' 
#' FDR functions for permutation-based estimators, including pi0 as well as FDR
#' confidence intervals.  The confidence intervals account for dependencies
#' between tests by the incorporation of an overdispersion parameter, which is
#' estimated from the permuted data.
#' 
#' \tabular{ll}{ Package: \tab fdrci\cr Type: \tab Package\cr Version: \tab
#' 2.1\cr Date: \tab 2018-02-21\cr License: \tab Artistic-2.0\cr LazyLoad: \tab
#' yes\cr } This method is designed to compute FDR when a permutation-based
#' approach has been utilized. The objective here is to identify a subset of
#' positive tests that have corresponding statistics with a more exteme
#' distribution than the permuted results, which are assumed to represent the
#' null. The significance of the subset is described in terms of the FDR and
#' uncertainty in the FDR estimate by computing a confidence interval. Say a
#' set of p-values(or simply a set of test statistics) were recorded for a set
#' of hypothesis tests, and data were permuted B times with test results
#' generated for each permutation. The function fdr_od() can be used to
#' estimate FDR and and a confidence interval along with pi0, the proportion of
#' true null hypotheses, given a selected significance threshold. The function
#' fdrTbl()uses fdr_od() to create a table of results over a sequence of
#' possible significance thresholds. Finally, the function FDRplot will plot
#' results from fdrTbl(), facilitating the selection of a final significance
#' threshold.
#' 
#' @name fdrci-package
#' @aliases fdrci-package fdrci
#' @docType package
#' @author Joshua Millstein
#' 
#' Maintainer: Joshua Millstein <millsteinjoshua@@gmail.com> Joshua Millstein
#' @references Millstein J, Volfson D. 2013. Computationally efficient
#' permutation-based confidence interval estimation for tail-area FDR.
#' Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
#' @keywords htest nonparametric
NULL



