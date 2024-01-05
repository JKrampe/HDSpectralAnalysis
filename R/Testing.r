# Coherence and Partial Coherence Testing ---------------------------------



#
#' The Cumulative distribution function of d independent chi2 distributions
#'
#' The Cumulative distribution function of d independent \eqn{chi^2}-distributed random variables
#'
#' @param t vector of value at which probabilities shall be computed
#' @param d Number of independent chi2 random variables
#'
#' @return probability at t
#' @export
#'
#' @examples
Gt2t=function(t,d)
{
  return(1-(pexp(t,rate=0.5))^d)
}

#' Helper function to determine threshold for test statistic
#'
#' This is a helper function to determine the threshold for the test statistic (TX) with asymptotic distribution (Gt).
#' Given q conducted tests, it gives the ratio between the expected number of positive test (under H0, i.e., distribution Gt) at threshold t
#' and the number positive test at threshold t
#'
#' @param t possible threshold t
#' @param TX Test statistic
#' @param Gt Asymptotic distribution of test statistic
#' @param q Number of test
#'
#' @return Given q conducted tests, it gives the ratio between the expected number of positive test (under H0, i.e., distribution Gt) at threshold t
#' and the number positive test at threshold t
#' @export
#'
#' @examples
t_inf_help=function(t,TX,Gt,q)
{
  return(Gt(t)*q/max(1,sum(TX[TX<Inf]>t)))
}


#' Computes threshold to FDR alpha
#'
#' @param TX Test statistic
#' @param Gt Asymptotic distribution of test statistic
#' @param q Number of test
#' @param ttilde Upper threshold to which point the asymptotic distribution of test statistic is trustworthy
#' @param lengt.out Number of threshold to try in interval (0,ttilde)
#' @param alpha desired FDR
#'
#' @return Derived threshold
#' @export
#'
#' @examples
t_inf=function(TX,Gt,q,ttilde,lengt.out,alpha)
{
  t_seq=seq(0,ttilde,length.out = lengt.out)[-1]
  return(c(t_seq,ttilde)[which.max(c(sapply(t_seq,function(t) t_inf_help(t,TX,Gt,q))<=alpha,T))])
}


#' Computes the empirical false discovery rate
#'
#' Computes the empirical false discovery rate for detecting (partial) coherences above the threshold Delta_0
#'
#' @param TX Test statistic
#' @param t_hat Chosen threshold
#' @param True_R True (partial) coherence
#' @param Delta_0 Delta in H0, i.e., (partial) coherence above that delta shall be identified
#'
#' @return empirical FDR
#' @export
#'
#' @examples
FDR_TX=function(TX,t_hat,True_R,Delta_0)
{
  return((sum(TX[True_R<=Delta_0]>t_hat))/max(1,(sum(TX>t_hat))))
}

#' Computes the empirical power
#'
#' Computes the empirical power for detecting (partial) coherences above the threshold Delta_0
#'
#' @param TX Test statistic
#' @param t_hat Chosen threshold
#' @param True_R True (partial) coherence
#' @param Delta_0 Delta in H0, i.e., (partial) coherence above that delta shall be identified
#'
#' @return empirical power
#' @export
#'
#' @examples
Emp_power=function(TX,t_hat,True_R,Delta_0)
{
  return((sum(TX[(True_R>Delta_0)]>t_hat))/(sum(True_R>Delta_0)))
}
