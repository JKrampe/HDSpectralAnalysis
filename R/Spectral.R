# Spectral Density Estimator ----------------------------------------------


#' S_comp
#' Given discrete Fourier transformed Z, computes \eqn{Z Z^*} and splits it into real and imaginary parts.
#'
#' @param Zt complex matrix of size n times p
#'
#' @return real matrix of size 2p times 2p
#'
#' @examples
S_comp=function(Zt)
{
  A1=var(Re(Zt))+var(Im(Zt))
  B1=cov(Im(Zt),Re(Zt))-cov(Re(Zt),Im(Zt))
  return(rbind(cbind(A1,B1),cbind(-B1,A1)))
}


#'  Adaptive rule to chose the lag-window
#'
#'  Adaptive rule to chose the lag-window, see Politis (2003): Adaptive bandwidth choice. This code is a modification of iosmooth::bwadap.ts by T. McMurry.
#'  here, the Frobenius-norm of the autocorrelation matrix are considered instead of the autocorrelation itself. The threshold is scaled by a factor p.
#'
#' @param X Multivariate time series
#' @param Kn Tuning parameter Kn discussed in Politis (2003). Roughly, the number of lags the autocorrelation function must stay below the threshold determined by c.thres
#' @param c.thres The bandwidth is chosen by looking for the first time the Frobenius norm of the autocorrelation function drops below \code{c.thres * p*sqrt(1/n)} and stays below that level for \code{Kn} lags.
#' @param acf.max Maximum number of autocorrelations to consider
#' @param even flag if always and even lag window should be returned
#' @param Kern used kernel
#'
#' @return list: M: suggested lag-window, M_corr_fac: Scaling factor for used kernel considered Jenkins kernel-bandwidth adjustment
#' @export
#'
#' @examples
adaptive_bandwidth_TLM=function(X,Kn=5,c.thres=2,acf.max=(dim(X)[1])^(1/2),even=T,Kern=Kernel)

{

  n=dim(X)[1]
  p=dim(X)[2]
  thresh <- c.thres * p*sqrt(log(n, 10)/n)
  ac <- apply(acf(Re(X), type = "correlation", plot = FALSE,
                  lag.max = min(acf.max,ceiling(n^(2/3)),floor(n/2)))$acf,1,norm,type="F")
  tmp=abs(ac)<thresh
  tmp2=tmp
  for(i in 1:(Kn-1))
  {
    tmp=c(tmp[-1],1)
    tmp2=tmp2+tmp
  }
  M=which.max(c(tmp2,Kn)==Kn)
  M1=ceiling(M/2)*2
  FourierFrequencies=0:(n-1)*2*pi/n
  #Jenkins kernel-bandwidth adjustment
  M_corr_fac=sum(Uniform(FourierFrequencies,M)^2)/(sum(Uniform(FourierFrequencies,M))^2)*M1/M/(sum(Kern(FourierFrequencies,M1)^2)/(sum(Kern(FourierFrequencies,M1))^2))
  M=M*M_corr_fac
  if(even==T)
    return(c(M=2*ceiling(M/2),corr_fac=M_corr_fac))
  return(c(M=ceiling(M),M_corr_fac))
}


#' Spectral density matrix inverse estimator
#'
#'  Discrete Fourier transformed Z are used to compute the spectral density matrix inverse f^{-1} at frequency Frequencies.
#'   A graphical lasso (glassofast) is used and tuning parameter lambda is chosen by (e)BIC
#'   Kernel smoothing is used with kernel Kernel and bandwith M
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param lambda tuning parameter for graphical lasso
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param Frequencies frequencies at which f^{-1} should be computed
#' @param THR flag if output should be thresholded at lambda
#' @param Check_pos_def flag if output is transformed to positive definiteness
#' @param M_corr Jenkins kernel-bandwidth adjustmentl. I.e., scaling to same variance as uniform kernel
#' @param gamma tuning parameter for ebic. Gamma=0 results in classical BIC
#' @param ... optional parameters for graphical lasso.
#'
#' @return list with array of f^{-1} of size p times p times length(Frequencies) and optimal (eBIC) lambdas.
#' @export
#'
#' @examples
spec_inv_freq=function(Z,lambda,Kernel,M,Frequencies,THR=FALSE,Check_pos_def=T,M_corr=1,gamma=0.5,...)
{

  # See equation (eq.fhat.inverse) for details.
  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n

  if(is.null(M))
  {
    M=ceiling(quantile(as.vector(adaptive_bandwidth_TLM(Z_to_X(Z))),probs = 0.8))
  }

  if(is.null(lambda))
  {
    lambda=sqrt(log(p)/n)*exp(seq(from = -15,to = 10,length.out =  max(p,n)))
  }
  erg=sapply(Frequencies,function(omega)
  {

    Weights=Kernel(omega-FourierFrequencies,M=M)
    S=S_comp((sqrt(Weights)*Z))
    if(Check_pos_def==T)
    {
      tmp1=eigen(S,symmetric = T)
      if(min(Re(tmp1$values))<1/n)
      {
        S=tmp1$vectors%*%diag(pmax(1/n,Re(tmp1$values)))%*%Conj(t(tmp1$vectors))
      }
    }

    tmp=glasso_BIC(S=S,lam=lambda,n=n/M*M_corr,gamma=gamma)

    lam=tmp$best_lam
    tmp=tmp$Omega
    erg_omg=((tmp[1:p,1:p]+tmp[1:p+p,1:p+p])+1i*(tmp[1:p,1:p+p]-tmp[1:p+p,1:p]))/2

    erg_omg=(erg_omg+Conj(t(erg_omg)))/2
    if(Check_pos_def==T)
    {
      tmp1=eigen(erg_omg,symmetric = T)
      if(min(Re(tmp1$values))<1/n)
      {
        erg_omg=tmp1$vectors%*%diag(pmax(1/n,Re(tmp1$values)))%*%Conj(t(tmp1$vectors))
      }
    }
    return(list(spec=erg_omg,lam=lam))
  })

  return(list(spec=simplify2array(erg[1,]),lambdas=simplify2array(erg[2,]),M=M))
}




#' Spectral density matrix inverse estimator
#'
#' Discrete Fourier transformed Z are used to compute the spectral density matrix inverse f^{-1} at frequency Frequencies.
#' A graphical lasso (glassofast) is computed for every tuning parameter lambda.
#' Kernel smoothing is used with kernel Kernel and bandwith M
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param lambda tuning parameter for graphical lasso
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param Frequencies frequencies at which f^{-1} should be computed
#'
#' @return array of f^{-1} of size p times p times length(Frequencies) times length (lambda)
#' @export
#'
#' @examples
spec_inv_freq_lam=function(Z,lambda,Kernel,M,Frequencies)
{
  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n

  if(is.null(M))
  {
    M=ceiling(quantile(as.vector(adaptive_bandwidth_TLM(Z_to_X(Z))),probs = 0.8))
  }

  if(is.null(lambda))
  {
    lambda=sqrt(log(p)/n)*exp(seq(from = -15,to = 10,length.out =  max(p,n)))
  }

  return(list(spec=sapply(Frequencies,function(omega)
  {
    Weights=Kernel(omega-FourierFrequencies,M=M)
    Weights=sqrt(Weights)
    Zt=((Weights)*Z)

    return(glasso_seq_lam(S=S_comp(Zt),lam=lambda,n=n/M)$Omega_seq)
  },simplify = "array"),M=M))
}





#' Spectral density matrix inverse estimator with prewhiteing
#'
#' An estimator for the spectral density inverse matrix. For this, first a sparse VAR is fitted to the data. Then on the residuals of the sparse VAR,
#' a nonparameteric spectral density estimator using graphical lasso with eBIC selection is fitted to compute the spectral density matrix inverse of the residuals.
#' Kernel smoothing is used with kernel Kernel and bandwith M.
#' Final output is the computed spectral density matrix inverse estimator. I.e, an estimator of \eqn{f^{-1}} of the multivariate time series X_t.
#'
#' @param X_t multivariate time sereis, matrix n times p
#' @param lambda tuning parameter for graphical lasso
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param M_corr Jenkins kernel-bandwidth adjustmentl. I.e., scaling to same variance as uniform kernel
#' @param VAR_order order of the used VAR model
#' @param Frequencies frequencies at which f^{-1} should be computed
#' @param crit used information criteria to select the tuning parameter lambda in the sparse VAR estimation
#' @param Check_pos_def flag if output is transformed to positive definiteness
#' @param lambda_VAR tuning parameter lambda in the sparse VAR estimation
#' @param gamma tuning parameter for ebic. Gamma=0 results in classical BIC
#' @param ... optional parameter pased to the sparse VAR estimation
#'
#'
#' @return list with f_inv: array of f^{-1} of size p times p times length(Frequencies); VAR_A: VAR model; f_eps_inv f^{-1} of the VAR residuals; M: bandwidth of f_eps_inv.
#' @export
#'
#' @examples
spec_inv_preWhiteningVAR_freq=function(X_t,lambda,Kernel,M,M_corr=1,VAR_order,Frequencies,crit="bic",Check_pos_def=T,lambda_VAR=lambda,gamma=0.5,...)
{
  n=dim(X_t)[1]
  p=dim(X_t)[2]
  if(is.null(p))
  {
    return(fhatcv2(lambda = Frequencies,X_t = X_t)^{-1})
  }
  if(VAR_order==0)
  {
    f_inv=spec_inv_freq(Z = X_to_Z(X_t),Frequencies = Frequencies,lambda = lambda,Kernel = Kernel,M = M,Check_pos_def=Check_pos_def)
    return(list(f_inv=f_inv$spec,M=f_inv$M,VAR_A=NULL))
  }
  m=length(Frequencies)
  VAR_est=VARLASSO(X_t,VAR_order=VAR_order,lambda = lambda_VAR,crit=crit,...)
  Z_eps=X_to_Z(VAR_est$eps)
  f_inv=spec_inv_freq(Z = Z_eps,Frequencies = Frequencies,lambda = lambda,Kernel = Kernel,M = M,Check_pos_def=Check_pos_def,M_corr=M_corr,gamma=gamma)
  Spec_eps=f_inv$spec

  A=A_to_array(VAR_est$A)
  tmp=sapply(1:VAR_order,function(i) matrix(i,p,p),simplify = "array")
  return(list(f_inv=sapply(1:m,function(omg)
    t(diag(p)-apply(A*exp(1i*tmp*Frequencies[omg]),1:2,sum))%*%Spec_eps[,,omg]%*%((diag(p)-apply(A*exp(-1i*tmp*Frequencies[omg]),1:2,sum))),simplify = "array"),
    VAR_A=A,f_eps_inv=Spec_eps,M=f_inv$M))
}







#' De-biased partial coherence with pre-filtering
#'  Computes for frequencies Frequencies the de-biased partial coherence (rho) for all u,v.
#'  A VAR-prefiltering (Filter) is used in the kernel estimation.
#'  spec should be an estimate of f^{-1} at frequencies Frequencies.
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param spec estimated f^{-1} at frequencies Frequencies of size p times p times length(Frequencies)
#' @param Filter Filter for the pre-filtering step. Filter input as a VAR filter as matrix of size p times (p*VARorder)
#' @param Frequencies Frequencies at which the de-biased partial coherence should be computed
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param DEMEAN flag if Z should be centered prior to kernel estimation
#'
#' @return rho_{u,v}(omega), omega in Frequencies, u,v=1,...,p  as array of p times p times length(Frequencies)
#' @export
#'
#' @examples
part_coh_de_Filter=function(Z,spec,Filter,Frequencies,Kernel,M,DEMEAN=FALSE)
{


  if(dim(spec)[3]!=length(Frequencies))
    stop("Spectral inv estimate and frequences do not match ")
  eps=VAR_Residuals(X_t = Z_to_X(Z),A = Filter)
  Z=X_to_Z(eps)
  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n
  A=A_to_array(Filter)
  tmp_A=sapply(1:dim(A)[3],function(i) matrix(i,p,p),simplify = "array")


  part_rho=sapply(1:length(Frequencies),function(omega_i)
  {
    Weights=Kernel(-FourierFrequencies+Frequencies[omega_i],M=M)
    Zt=Z*Weights
    if(DEMEAN)
      Zt=t(t(Zt)-apply(Zt,2,mean))
    Filter_w=solve((diag(p)-apply(A*exp(-1i*tmp_A*Frequencies[omega_i]),1:2,sum)))
    tmp=Filter_w%*%t(Zt)%*%Conj(Z)%*%t(Conj(Filter_w))%*%spec[,,omega_i]
    tmp1=spec[,,omega_i]%*%tmp
    spec_diag=Re(diag(spec[,,omega_i]))
    Bias_cor=sapply(1:p,function(u) sapply(1:p,function(v) tmp1[v,u]-tmp1[v,v]*spec[v,u,omega_i]/spec[v,v,omega_i]))
    DN=sapply(1:p,function(u) sapply(1:p,function(v) tmp[u,u]*spec[v,v,omega_i]-tmp[u,v]*spec[v,u,omega_i]))


    tmp2=-(diag(1/spec_diag)%*%t(spec[,,omega_i]))+(Bias_cor/DN)
    diag(tmp2)=1
    beta_de_scaled=tmp2*outer(sqrt(spec_diag),Y = sqrt(1/spec_diag),FUN = "*")
    return(0.5*Conj(beta_de_scaled+t(Conj(beta_de_scaled))))
  },simplify = "array")

  return(part_rho)
}


#' De-biased partial coherence
#'  Computes for frequencies Frequencies the de-biased partial coherence (rho) for all u,v.
#'  spec should be an estimate of f^{-1} at frequencies Frequencies.
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param spec estimated f^{-1} at frequencies Frequencies of size p times p times length(Frequencies)
#' @param Frequencies Frequencies at which the de-biased partial coherence should be computed
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param DEMEAN flag if Z should be centered prior to kernel estimation
#'
#' @return rho_{u,v}(omega), omega in Frequencies, u,v=1,...,p  as array of p times p times length(Frequencies)
#' @export
#'
#' @examples
part_coh_de_red=function(Z,spec,Frequencies,Kernel,M,DEMEAN=FALSE)
{
  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n

  part_rho=sapply(1:length(Frequencies),function(omega_i)
  {
    Weights=Kernel(-FourierFrequencies+Frequencies[omega_i],M=M)

    Zt=Z*Weights
    if(DEMEAN)
      Zt=t(t(Zt)-apply(Zt,2,mean))

    tmp=t(Zt)%*%Conj(Z)%*%spec[,,omega_i]
    tmp1=spec[,,omega_i]%*%tmp
    spec_diag=Re(diag(spec[,,omega_i]))
    Bias_cor=sapply(1:p,function(u) sapply(1:p,function(v) tmp1[v,u]-tmp1[v,v]*spec[v,u,omega_i]/spec[v,v,omega_i]))
    DN=sapply(1:p,function(u) sapply(1:p,function(v) tmp[u,u]*spec[v,v,omega_i]-tmp[u,v]*spec[v,u,omega_i]))


    tmp2=-(diag(1/spec_diag)%*%t(spec[,,omega_i]))+(Bias_cor/DN)
    diag(tmp2)=1
    beta_de_scaled=tmp2*outer(sqrt(spec_diag),Y = sqrt(1/spec_diag),FUN = "*")
    return(0.5*Conj(beta_de_scaled+t(Conj(beta_de_scaled))))
  },simplify = "array")
}




#' Chi^2 version of de-biased partial coherence
#'
#' Computes the Chi^2 version of de-biased partial coherence. For this, first the variance matrix of the real and imaginary part is computed.
#' Then, the inverse matrix of this is used as a transformtion on the real and imaginary parts.
#'
#' @param x vectorized array of partial coherences of size p times p times d
#' @param y vectorized array of test statistics of size p times p times d
#'
#' @return Chi^2 version of de-biased partial coherence
#'
#' @examples
var_part_coherence_inv=function(x,y)
{
  sig1=(1-Re(x)^2)
  sig2=(1-Im(x)^2)
  sig3=(Re(x)*Im(x))
  return(2*(Re(y)^2*sig2+sig3*Re(y)*Im(y)+sig1*Im(y)^2)/(pmax(10E-5,((1-abs(x)^2)^2))))
}



#' Test statistic based on the maximal de-biased partial coherence with pre-filtering
#'
#' Computes the test statistic based on the maximum partial coherences over the frequencies (FourierFrequencies_h).
#' The partial coherence is computed with the pre-filter de-biased partial coherence estimator.
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param spec estimated f^{-1} at frequencies Frequencies of size p times p times length(Frequencies)
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param DEMEAN flag if Z should be centered prior to kernel estimation
#' @param Delta_0 specify Delta 0 is test is not performed to Delta_0=0.
#' @param FourierFrequencies_h Fourier Frequencies at which the de-biased partial coherence should be computed
#' @param Filter Filter for the pre-filtering step. Filter input as a VAR filter as matrix of size p times (p*VARorder)

#'
#' @return list: TX_Z: test statistic in raw form; TX_chi: test statistic in chi^2 form; part_rho: de-biased partial coherence; C_K: kernel variance constant; M: used bandwidth
#' @export
#'
#' @examples
max_de_part_coh_H_0_pre_Filt=function(Z,spec,Kernel,M,DEMEAN=TRUE,Delta_0,FourierFrequencies_h,Filter)
{

  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n

  nh=(n/M)


  C_K=sum(Kernel(FourierFrequencies,M)^2)/(sum(Kernel(FourierFrequencies,M))^2)*nh


  erg_part=part_coh_de_Filter(Z = Z,spec=spec,Frequencies = FourierFrequencies_h,Kernel = Kernel,M = M,DEMEAN = DEMEAN,Filter = Filter)

  if(!is.null(Delta_0))
  {
    TX_tmp=as.numeric(abs(erg_part)>Delta_0)*(erg_part-Delta_0*exp(1i*Arg(erg_part)))
  }else{
    TX_tmp=erg_part
  }

  spec_rho=f_inv_to_part_rho(spec)
  TX_chi=array(var_part_coherence_inv(x = as.vector(spec_rho),y=as.vector(TX_tmp)),dim=dim(TX_tmp))/C_K
  var_freq=C_K*(1-abs(spec_rho)^2)*(1-(1/2+1/2*as.numeric(FourierFrequencies_h==0))*abs(spec_rho)^2)
  TX_Z=abs(TX_tmp)^2/var_freq


  return(list(TX_Z=apply(TX_Z,1:2,max),TX_chi=apply(TX_chi,1:2,max),part_rho=erg_part,C_K=C_K,M=M))
}


#' Test statistic based on the maximal coherence with pre-filtering
#'
#' Computes the test statistic based on the maximum coherences over the frequencies (FourierFrequencies_h).
#' The coherence is computed with the pre-filtered lag-window coherence estimator.
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param DEMEAN flag if Z should be centered prior to kernel estimation
#' @param Delta_0 specify Delta 0 is test is not performed to Delta_0=0.
#' @param FourierFrequencies_h Fourier Frequencies at which the de-biased partial coherence should be computed
#' @param Filter Filter for the pre-filtering step. Filter input as a VAR filter as matrix of size p times (p*VARorder)
#'
#' @return list: TX_Z: test statistic in raw form; TX_chi: test statistic in chi^2 form; coh_rho: partial coherence; C_K: kernel variance constant; M: used bandwidth
#' @export
#'
#' @examples
max_coh_H_0_pre_Filt=function(Z,Kernel,M,DEMEAN=TRUE,Delta_0,FourierFrequencies_h,Filter)
{
  eps=VAR_Residuals(X_t = Z_to_X(Z),A = Filter)
  Z=X_to_Z(eps)
  n=dim(Z)[1]
  p=dim(Z)[2]

  nh=(n/M)
  FourierFrequencies=0:(n-1)*2*pi/n
  C_K=sum(Kernel(FourierFrequencies,M)^2)/(sum(Kernel(FourierFrequencies,M))^2)*nh
  A=A_to_array(Filter)
  tmp_A=sapply(1:dim(A)[3],function(i) matrix(i,p,p),simplify = "array")


  erg_coh=sapply(1:length(FourierFrequencies_h),function(omega_i)
  {
    Weights=Kernel(-FourierFrequencies+FourierFrequencies_h[omega_i],M=M)
    Zt=Z*Weights
    if(DEMEAN)
      Zt=t(t(Zt)-apply(Zt,2,mean))
    Filter_w=solve((diag(p)-apply(A*exp(-1i*tmp_A*FourierFrequencies_h[omega_i]),1:2,sum)))
    tmp=Filter_w%*%t(Zt)%*%Conj(Z)%*%t(Conj(Filter_w))
    d_tmp=sqrt(Re(diag(tmp)))
    return(t(t(tmp/d_tmp)/d_tmp))
  },simplify = "array")


  if(!is.null(Delta_0))
  {
    TX_tmp=as.numeric(abs(erg_coh)>Delta_0)*(erg_coh-Delta_0*exp(1i*Arg(erg_coh)))
  }else{
    TX_tmp=erg_coh
  }

  TX_chi=array(var_part_coherence_inv(x = as.vector(erg_coh),y=as.vector(TX_tmp)),dim=dim(TX_tmp))/C_K
  var_freq=C_K*(1-Delta_0^2)*(1-(1/2+1/2*as.numeric(FourierFrequencies_h==0))*Delta_0^2)
  TX_Z=abs(TX_tmp)^2/var_freq


  return(list(TX_Z=apply(TX_Z,1:2,max),TX_chi=apply(TX_chi,1:2,max),coh_rho=erg_coh,C_K=C_K,M=M))
}





#' Test statistic based on the maximal de-biased partial coherence
#'
#' Computes the test statistic based on the maximum partial coherences over the frequencies (FourierFrequencies_h).
#' The partial coherence is computed with the de-biased partial coherence estimator.
#'
#'
#' @param Z Discrete Fourier transformed as matrix n (sample size) times p (dimension)
#' @param spec estimated f^{-1} at frequencies Frequencies of size p times p times length(Frequencies)
#' @param Kernel used kernel in frequency domain as a function with input frequencies and lag window
#' @param M lag window for kernel
#' @param DEMEAN flag if Z should be centered prior to kernel estimation
#' @param Delta_0 specify Delta 0 is test is not performed to Delta_0=0.
#' @param FourierFrequencies_h Fourier Frequencies at which the de-biased partial coherence should be computed

#'
#' @return list: TX_Z: test statistic in raw form; TX_chi: test statistic in chi^2 form; part_rho: de-biased partial coherence; C_K: kernel variance constant; M: used bandwidth
#' @export
max_de_part_coh_H_0=function(Z,spec,Kernel,M,DEMEAN=TRUE,Delta_0,FourierFrequencies_h)
{
  # Computes the maximum Test statistic (eq.test) for u,v=1,...,p
  # to threshold Delta_0
  # spec should be an estimate of f^{-1} at frequencies Frequencies
  ## Teststatistic for all u,v=1,...,p
  n=dim(Z)[1]
  p=dim(Z)[2]
  FourierFrequencies=0:(n-1)*2*pi/n
  nh=(n/M)
  C_K=sum(Kernel(FourierFrequencies,M)^2)/(sum(Kernel(FourierFrequencies,M))^2)*nh

  erg_part=part_coh_de_red(Z = Z,spec=spec,Frequencies = FourierFrequencies_h,Kernel = Kernel,M = M,DEMEAN = DEMEAN)

  if(!is.null(Delta_0))
  {
    TX_tmp=as.numeric(abs(erg_part)>Delta_0)*(erg_part-Delta_0*exp(1i*Arg(erg_part)))
  }else{
    TX_tmp=erg_part
  }

  spec_rho=f_inv_to_part_rho(spec)
  TX_chi=array(var_part_coherence_inv(x = as.vector(spec_rho),y=as.vector(TX_tmp)),dim=dim(TX_tmp))/C_K
  var_freq=C_K*(1-abs(spec_rho)^2)*(1-(1/2+1/2*as.numeric(FourierFrequencies_h==0))*abs(spec_rho)^2)
  TX_Z=abs(TX_tmp)^2/var_freq


  return(list(TX_Z=apply(TX_Z,1:2,max),TX_chi=apply(TX_chi,1:2,max),part_rho=erg_part,C_K=C_K,M=M))
}






#' Spectral density matrix inverse at frequency Omega for a VAR(d) model
#'
#' Compute the inverse spectral density at frequences Omega for a VAR(d) model with parameter A and variance matrix of the residuals Sigma.
#' Either Sigma or Sigma inverse is required.
#'
#' @param omega vector of frequencies
#' @param A VAR model matrix as matrix p times p*VARorder
#' @param Sigma_inv Sigma inverse matrix as p times p matrix
#' @param Sigma Sigma  matrix as p times p matrix
#'
#' @return spectral density inverse matrix as array of size p times p times length(Omega)
#' @export
#'
#' @examples
Spec_inv_VAR=function(omega,A,Sigma_inv=NULL,Sigma=NULL)
{

  #Dimension p
  #A is given a a p x p*d
  p=dim(A)[1]
  d=dim(A)[2]/p
  if(is.null(Sigma_inv))
    Sigma_inv=solve(Sigma)
  A=A_to_array(A)
  tmp=sapply(1:d,function(i) matrix(i,p,p),simplify = "array")
  return(sapply(omega,function(omg)
    2*pi*t(diag(p)-apply(A*exp(1i*tmp*omg),1:2,sum))%*%Sigma_inv%*%((diag(p)-apply(A*exp(-1i*tmp*omg),1:2,sum))),simplify = "array"))
}

##' Spectral density matrix inverse at frequency Omega for a VMA(d) model
#'
#' Compute the inverse spectral density at frequences Omega for a VMA(d) model with parameter A and variance matrix of the residuals Sigma.
#'
#' @param omega vector of frequencies
#' @param A VMA model matrix as matrix p times p*MA order
#' @param Sigma Sigma  matrix as p times p matrix
#'
#' @return spectral density inverse matrix as array of size p times p times length(Omega)
#' @export
#'
#' @examples
f_spec_MA=function(A,Sigma,omega_list)
{
  p=dim(A)[1]
  d=dim(A)[2]/p
  A_array=A_to_array(A)
  return(sapply(omega_list,function(omega)
  {
    A_1=diag(d)+apply(sapply(1:d,function(i) A_array[,,i]*exp(-1i*i*omega),simplify = "array"),1:2,sum)

    return(1/(2*pi)*A_1%*%Sigma%*%t(Conj(A_1)))
  },simplify = "array"))
}
##' Spectral density matrix inverse at frequency Omega for a VARMA(d,q) model
#'
#' Compute the inverse spectral density at frequences Omega for a VARMA(d,q) model with parameter A (VAR), B (VMA) and variance matrix of the residuals Sigma.
#'
#' @param A VAR model matrix as matrix p times p*VAR order
#' @param B VMA model matrix as matrix p times p*VMA order
#' @param Sigma Sigma  matrix as p times p matrix
#' @param omega_list vector of frequencies
#'
#' @return spectral density inverse matrix as array of size p times p times length(omega_list)
#' @export
#'
#' @examples
f_spec_inv_ARMA=function(A,B,Sigma,omega_list)
{
  p=dim(A)[1]
  d=dim(A)[2]/p
  d2=dim(B)[2]/p
  A_array=A_to_array(A)
  B_array=A_to_array(B)
  Sigmainv=solve(Sigma)
  tmp=sapply(1:d,function(i) matrix(i,p,p),simplify = "array")
  return(sapply(omega_list,function(omega)
  {
    A_1=solve(diag(p)+apply(sapply(1:d2,function(i) B_array[,,i]*exp(1i*i*omega),simplify = "array"),1:2,sum))

    return(2*pi*t(diag(p)-apply(A_array*exp(1i*tmp*omega),1:2,sum))%*%t(A_1)%*%Sigmainv
           %*%Conj(A_1)%*%((diag(p)-apply(A_array*exp(-1i*tmp*omega),1:2,sum))))
  },simplify = "array"))
}


