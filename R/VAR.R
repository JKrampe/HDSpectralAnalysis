# Sparse VAR Models -------------------------------------------------------

#' Estimates a sparse VAR model
#'
#' Estimates a sparse VAR model using a row-wise lasso (glmnet). The tuning parameter in the lasso estimation can be chosen by various information criteria, see also see also \code{\link{Regression_row}}
#'
#'
#' @param X_t multivariate times series as matrix n times p
#' @param lambda vector of tuning paramter lambda
#' @param THR flag to specify if thresholded (with lambda as threshold) VAR model should be returned.
#' @param SCALE flag to specify if X should be normalized prior to regularized regression
#' @param ADAPTIVE flag to specify if a second adaptive step should be performed
#' @param crit flag to specify which information criteria should be used
#' @param tau if a second adptive is used, tau specifies the following penalty.factor=1/abs(abs(first.step.coef)+1/sqrt(n))^tau
#' @param first.penalty.factor penalty factors for the first step. NULL leads to uninformative first step factors.
#' @param VAR_order Order of the VAR model
#' @param ... additional parameters passed to glmnet
#' @param mbic2_c additinal parameter for some information criteria
#' @param Variablen vector of indices for which the row-wise VAR model should be computed if interest is only in a subset.
#'
#' @return list: A: sparse VAR model as matrix p times p*VARorder; eps: VAR residuals as matrix (n-VAR_order)times p; lambda: selected lambda by IC as vector of length p.
#' @export
#'
#' @examples
VARLASSO=function(X_t,VAR_order,Variablen=NULL,lambda,THR=FALSE,SCALE=TRUE,ADAPTIVE=TRUE,tau=1,mbic2_c=4,crit = c("bic", "aic", "aicc", "hqc","mbic2","eric","bic_m"),first.penalty.factor=NULL,...)
{
  lambda=sort(lambda,decreasing = TRUE)
  n = dim(X_t)[1]
  p = dim(X_t)[2]
  if(VAR_order==0)
    return(list(A=NULL,eps=X_t,lambda=NULL))
  Xtd = matrix(0, n - VAR_order, VAR_order * p)
  for (j in 1:VAR_order)
    Xtd[, (j - 1) * p + 1:p] = as.matrix(X_t[(VAR_order + 1 - j):(n - j), ])
  Y=X_t[(VAR_order + 1):n,]

  A_erg=matrix(0,p,p*VAR_order)

  lambda_sel=double(p)
  Residual=matrix(0,n-VAR_order,p)

  if(is.null(Variablen))
    Variablen=1:p

  for(j in Variablen)
  {
    MODEL_SELECT=Regression_row(Y = Y[,j],X = Xtd,lambda = lambda,THR = THR,SCALE = SCALE,ADAPTIVE = ADAPTIVE,crit = crit,tau = tau,mbic2_c = mbic2_c,first.penalty.factor=first.penalty.factor,...)

    A_erg[j,]=MODEL_SELECT$coef
    Residual[,j]=MODEL_SELECT$eps
    lambda_sel[j]=MODEL_SELECT$lam
  }

  return(list(A=A_erg,eps=Residual,lambda=lambda_sel))
}



#' Computes the Residuals of a VAR Model
#' Computes the Residuals of a VAR Model
#'
#' @param X_t Time Series X_t (n x p)
#' @param A VAR model parameter matrix A  as matrix p x (VAR_order * p)
#'
#' @return Residuals as matrix (n-VAR_order)times p
#' @export
#'
#' @examples
VAR_Residuals=function(X_t,A)
{
  n=dim(X_t)[1]
  p=dim(X_t)[2]
  VAR_order=dim(A)[2]/p
  Xtd = matrix(0, n - VAR_order, VAR_order * p)
  for (j in 1:VAR_order)
    Xtd[, (j - 1) * p + 1:p] = as.matrix(X_t[(VAR_order + 1 - j):(n - j), ])
  return(X_t[(VAR_order + 1):n,]-Xtd%*%t(A))
}



#' Simulates a VARMA model
#'
#' Generates n realization of a p-dimensional VARMA(VAR order,VMA order) model.
#'
#' @param A VAR model matrix as matrix p times p*VAR order
#' @param B VMA model matrix as matrix p times p*VMA order
#' @param Sigma Sigma  matrix as p times p matrix
#' @param n number of samples to generate
#' @param n2 burn-in period
#'
#' @return multivariate time series X_t as matrix n times p
#' @export
#'
#' @examples
simulateVARMAs=function(A,B, Sigma, n, n2 = 100)
{
  p = dim(A)[1]
  d = max(0,dim(A)[2] / p)
  d2 = max(0,dim(B)[2] / p)
  eps_t = mvtnorm::rmvnorm(n + n2, sigma = Sigma)
  Z_t=eps_t
  if(d2>=1)
    for (j in (d2 + 1):(n + n2))
      Z_t[j, ] =eps_t[j, ]+B%*% as.vector(sapply(1:d2, function(i)
        return(eps_t[j - i, ])))
  X_t=Z_t
  if(d>=1)
    for (j in max((d + 1),(d2 + 1)):(n + n2))
      X_t[j, ] =Z_t[j, ]+ A %*% as.vector(sapply(1:d, function(i)
        return(X_t[j - i, ])))
  return(X_t[-(1:n2), ])
}
