#' Given Regression Y=MODEL(i) X , computes several information criteria  and selects the best MODEL i
#'
#'   Given Regression Y=MODEL(i) X , computes several information criteria  and selects the best MODEL i.
#'   Information criteria: (n sample size, p dimension, n_var number of non-zero variables)
#'    bic = n * log(mse) + nvar * log(n)
#'    mbic2 = n * log(mse) + nvar * (log(n)+log(p/mbic2_c))-2*log(gamma(nvar+1))
#'    aic = n * log(mse) + 2 * nvar
#'    aicc = aic + (2 * nvar * (nvar + 1))/(n - nvar - 1)
#'    hqc = n * log(mse) + 2 * nvar * log(log(n))
#'    eric=n * log(mse) +nvar * log(sqrt(n*mse)/lambda)*ifelse(n>p,1,0.5)
#'    bic_m = n * log(mse) + nvar * log(n)*mbic2_c
#'
#' @param Y vector y of length n
#' @param X matrix X of size n times p
#' @param MODEL model parameter Beta ( Y=BETA X ) as array p times length(lambda)
#' @param crit flag to specify which information criteria should be used
#' @param lambda vector of tuning paramter lambda
#' @param mbic2_c additinal parameter for some information criteria
#' @return list: coef: Best model according to used IC; eps: Residuals of best model; lam: optimal lambda accoirding to used IC; crit: Minimal value of used IC.
#' @export
#'
#' @examples
Regression_Row_IC=function(Y,X,MODEL,crit = c("bic", "aic", "aicc", "hqc","mbic2","eric","bic_m"),lambda,mbic2_c=4)
{
  n = length(Y)
  nlambda=length(lambda)
  crit = match.arg(crit)
  df = apply(round(abs(MODEL),6)!=0,2,sum)
  p=dim(X)[2]
  nu_eric=ifelse(n>p,1,0.5)
  residuals = Y-X%*%MODEL
  mse=(apply(residuals,2,var)+as.numeric(df>(n/2))*(sqrt(1/n)*as.vector(var(Y))))

  sse = apply(residuals^2,2,sum)
  nvar = df
  bic = n * log(mse) + nvar * log(n)
  mbic2 = n * log(mse) + nvar * (log(n)+log(p/mbic2_c))-2*log(gamma(nvar+1))
  aic = n * log(mse) + 2 * nvar
  aicc = aic + (2 * nvar * (nvar + 1))/(n - nvar - 1)
  hqc = n * log(mse) + 2 * nvar * log(log(n))
  eric=n * log(mse) +nvar * log(sqrt(n*mse)/lambda)*nu_eric
  bic_m = n * log(mse) + nvar * log(n)*mbic2_c

  crit = switch(crit, bic = bic, aic = aic, aicc = aicc, hqc = hqc,mbic2=mbic2,eric=eric,bic_m=bic_m)

  selected=which.max(crit == min(crit))
  best.model=MODEL[,selected]
  best.residuals=residuals[,selected]

  return(list(coef=best.model,eps=best.residuals,lam=lambda[selected],crit=min(crit)))
}





#' Regularized linear regression Y ~ X
#'
#'  Computes the  regression Y=BETA X with l1 regularization. Uses LASSO (glmnet), a second adaptive step (penalty is chosen adaptively based on the first step (tau))
#'  Penalty parameter is chosen by information criteria (crit), see also \code{\link{Regression_Row_IC}}.
#'
#'  Information criteria: (n sample size, p dimension, n_var number of non-zero variables)
#'    bic = n * log(mse) + nvar * log(n)
#'    mbic2 = n * log(mse) + nvar * (log(n)+log(p/mbic2_c))-2*log(gamma(nvar+1))
#'    aic = n * log(mse) + 2 * nvar
#'    aicc = aic + (2 * nvar * (nvar + 1))/(n - nvar - 1)
#'    hqc = n * log(mse) + 2 * nvar * log(log(n))
#'    eric=n * log(mse) +nvar * log(sqrt(n*mse)/lambda)*ifelse(n>p,1,0.5)
#'    bic_m = n * log(mse) + nvar * log(n)*mbic2_c
#'
#' @param Y vector y of length n
#' @param X matrix X of size n times p
#' @param lambda vector of tuning paramter lambda
#' @param THR flag to specify if thresholded (with lambda as threshold) beta should be returned.
#' @param SCALE flag to specify if X should be normalized prior to regularized regression
#' @param ADAPTIVE flag to specify if a second adaptive step should be performed
#' @param crit flag to specify which information criteria should be used
#' @param mbic2_c additinal parameter for some information criteria
#' @param tau if a second adptive is used, tau specifies the following penalty.factor=1/abs(abs(first.step.coef)+1/sqrt(n))^tau
#' @param first.penalty.factor penalty factors for the first step. NULL leads to uninformative first step factors.
#' @param ... additional parameters
#'
#' @return list: coef: Best model according to used IC; eps: Residuals of best model; lam: optimal lambda accoirding to used IC; crit: Minimal value of used IC.
#' @export
#'
#' @examples
Regression_row=function(Y,X,lambda,THR=FALSE,SCALE=TRUE,ADAPTIVE=TRUE,crit = c("bic", "aic", "aicc", "hqc","mbic2","eric","bic_m"),mbic2_c=4,tau=1,first.penalty.factor=NULL,...)
{


  n = length(Y)
  p=dim(X)[2]
  nlambda=length(lambda)
  if(is.null(first.penalty.factor))
    first.penalty.factor=rep(1,p)
  crit = match.arg(crit)
  if(SCALE)
  {
    C=apply(X,2,sd)
    X=scale(X,center = FALSE,scale = TRUE)
  }else{
    C=double(p)+1
  }

  tmp=coef(glmnet::glmnet(x = X,y = Y,family = "gaussian",lambda = lambda,intercept = TRUE,standardize = FALSE,penalty.factor = first.penalty.factor,...))[-1,]

  tmp=(as.matrix(tmp)/C)
  MODEL_SELECT=Regression_Row_IC(Y = Y,X = t(t(X)*C),tmp,crit,lambda,mbic2_c = mbic2_c)
  if(ADAPTIVE)
  {
    first.step.coef = MODEL_SELECT$coef*C
    penalty.factor=1/abs(abs(first.step.coef)+1/sqrt(n))^tau
    tmp=coef(glmnet::glmnet(x = X,y = Y,family = "gaussian",lambda = lambda,intercept = TRUE,standardize = FALSE,penalty.factor = penalty.factor,...))[-1,]
    tmp=(as.matrix(tmp)/C)
    MODEL_SELECT=Regression_Row_IC(Y = Y,X = t(t(X)*C),tmp,crit,lambda,mbic2_c = mbic2_c)
  }
  return(MODEL_SELECT)

}




#' Graphical lasso with BIC
#'
#' This is a wrapper for the glassoFast function to compute the graphical lasso.
#' Additionally, the tuning parameter of the regularization parameter lambda is chosen by BIC (taking the following considerations into account: Wang, H., Li, B., and Leng, C. (2009). Shrinkage tuning parameter selection with a diverging number of parameters)
#' The (real-valued) input matrix S is standardized prior to computation.
#'
#' @param S Input matrix, i.e., empirical covariance matrix or (transformed) lag-window spectral density matrix
#' @param n (Effective) sample size (important for BIC)
#' @param lam regularization parameter lambda
#' @param maxit maximal number of iterations
#' @param tol Tolerance in numerical optimization
#' @param gamma Gamma-Parameter in BIC selection with considerations of Wang. et. al (2009)
#'
#' @return list: Lambdas: used lambda seqeuence, maxIt: maximal number of iterations,Omega: Estimated precision matrix or f^{-1},Sigma: Estimated Covariance matrix of f,BICs: BIC values, best_lam: BIC-optimal lambda
#' @export
#'
#' @examples
glasso_BIC=function(S,n,lam,maxit = 10000,tol = 1e-04,gamma=4)
{
  lam = sort(lam,decreasing = F)
  p=dim(S)[1]
  sds=sqrt(diag(S))
  S=diag(1/sds)%*%S%*%diag(1/sds)
  init = S
  initOmega = diag(ncol(S))
  BICs=rep(Inf,length(lam))
  best_S=NULL
  best_Omega=NULL
  best_lam=0
  for (i in 1:length(lam)) {
    lam_ = lam[i]

    GLASSO = glassoFast(S =S, rho = lam_, thr = tol, maxIt=maxit,start = "warm", w.init = init, wi.init = initOmega,
                        trace = FALSE)
    init = GLASSO$w
    initOmega = GLASSO$wi
    BICs[i]=min(Inf,(n) * 2*(sum(diag(GLASSO$wi%*%S)) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])+sum(GLASSO$wi[upper.tri(GLASSO$wi)] !=0) *( log(n)+gamma*log((p*(p-1)/2))),na.rm = T)
    if(BICs[i]==min(BICs))
    {
      best_S=GLASSO$w
      best_Omega=GLASSO$wi
      best_lam=lam_
    }
  }

  return(list(Lambdas = lam, maxIt=maxit, Omega = diag(1/sds)%*%best_Omega%*%diag(1/sds), Sigma = diag(sds)%*%best_S%*%diag(sds),BICs=BICs,best_lam=best_lam))
}

#' Graphical lasso
#'
#' This is a wrapper for the glassoFast function to compute the graphical lasso over a sequence of regularization parameters lambda.
#' The (real-valued) input matrix S is standardized prior to computation.
#'
#' @param S Input matrix, i.e., empirical covariance matrix or (transformed) lag-window spectral density matrix
#' @param n (Effective) sample size (important for BIC)
#' @param lam regularization parameter lambda
#' @param maxit maximal number of iterations
#' @param tol Tolerance in numerical optimization
#'
#' @return list: Lambdas: used lambda seqeuence, Omega_seq: Sequences (length(lam)) of estimated precision matrix or f^{-1}
#' @export
#' @export
#'
#' @examples
glasso_seq_lam=function(S,n,lam,maxit = 10000,tol = 1e-04)
{
  lam = sort(lam)
  sds=sqrt(diag(S))
  S=diag(1/sds)%*%S%*%diag(1/sds)
  init = S
  initOmega = diag(ncol(S))
  p=dim(S)[1]/2
  Omega_seq=array(0,c(p,p,length(lam)))
  for (i in 1:length(lam)) {
    lam_ = lam[i]

    GLASSO = glassoFast(S =S, rho = lam_, thr = tol, maxIt=maxit,start = "warm", w.init = init, wi.init = initOmega,
                        trace = FALSE)
    init = GLASSO$w
    initOmega = GLASSO$wi
    tmp= diag(1/sds)%*%GLASSO$wi%*%diag(1/sds)
    erg_omg=((tmp[1:p,1:p]+tmp[1:p+p,1:p+p])+1i*(tmp[1:p,1:p+p]-tmp[1:p+p,1:p]))/2
    erg_omg=(erg_omg+Conj(t(erg_omg)))/2
    Omega_seq[,,i]=erg_omg
  }

  return(list(Lambdas = lam, Omega_seq=Omega_seq))
}
