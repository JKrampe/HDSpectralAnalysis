#' Plot a matrix
#'
#' This functions plots a matrix A, i.e, x- axis are the column indices, y-axis are the row indices.
#' Color1 is used for positive, Color2 for negative entries. Zero entries are left blank. Entries closer to zero are plotted more transparent.
#'
#'
#' @param A matrix to plot
#' @param color1 Color for positive entries
#' @param color2 Color for negative entries
#' @param range_A Value range to consider
#' @param xlab Name for x-axis
#' @param ylab name for y-axis
#' @param pch pch-number used for plotted points (15 default)
#' @param alpha_func Function used for color-transparency
#' @param ... Additional parameters for the plot function
#'
#' @return No return value
#' @export
#'
#' @examples
plotmat = function(A,
                   color1 = "red",
                   color2 = "blue",
                   range_A = NULL,xlab = "column",ylab = "row",pch=15,alpha_func=function(x,y) scales::alpha(x,y),
                   ...)
{
  if (is.null(range_A))
    range_A = range(A)
  p1 = dim(A)[1]
  p2 = dim(A)[2]
  plot(
    c(1,p1),c(1,p2),
    col = "white",
    ylim = c(p1, 1),
    xlim = c(1, p2),
    xlab = xlab,
    ylab = ylab,
    ...
  )
  for (i in 1:p1)
    for (j in 1:p2)
      if (A[i, j] != 0)
      {
        if (A[i, j] > 0)
        {
          points(j,
                 i,
                 pch = pch,
                 col = alpha_func(color1, (A[i, j]/(range_A[2]))))
        } else
        {
          points(j,
                 i,
                 pch = pch,
                 col = alpha_func(color2, (A[i, j]/(range_A[1]))))
        }
      }
}



#' Computes the (partial) coherence network
#'
#' This function computes for a given time series X_t its coherence and partial coherence network. The objective is to identify all (partial) coherences which have values above Delta_0 over the frequency band (Frequencyband).
#' For this, first the spectral density matrix inverse is estimated using a pre-whitten (sparse VAR) lag-window estimator (with kernel Kernel) with graphical lasso regularization. Then, the de-biased partial coherence (with kernel Kernel) is computed and the test statistics
#' are threshold automatically to obtain (asymptotically) the FDR FDR_alpha.
#'
#'
#'
#' @param X_t time series
#' @param Frequencyband lower and upper value of the frequency band (subset of (0,pi)) to consider
#' @param FDR_alpha desired FDR
#' @param Delta_0  Delta in H0, i.e., (partial) coherence above that delta shall be identified. Delta0=0 means all non-zero (partial) coherences
#' @param Kernel Used kernel
#' @param r r is determined by the decay behavior of the Fourier coefficients of the kernel used; see Assumption 2 in Section 6. r=NULL for uniform kernel, r=2 for modified Bartlett
#' @param lambda lambda sequence for regularization
#'
#' @return list: Frequencyband_Fourier: Used Fourier frequencies; Partial_Coherence: Partial coherence network, Coherence: Coherence network
#' @export
#'
#' @examples
PC_Network=function(X_t,Frequencyband=c(0,pi),FDR_alpha=0.05,Delta_0=0,Kernel=Modified_Bartlett,r=2,lambda=NULL,VAR_order=ceiling(log(dim(X_t)[1],10)),bandwidth=NULL,mbic2_c=log(dim(X_t)[2]),gamma=log(log(dim(X_t)[2])))
{
  X_t=scale(X_t,center=T,scale=F)
  n=dim(X_t)[1]
  p=dim(X_t)[2]
  Z=X_to_Z(X_t)



  q=(p^2-p)/2
  q_ind=sapply(1:p,function(i) sapply(1:p,function(j) i>j),simplify = "matrix")

  if(is.null(lambda))
  {
    lambda=sqrt(log(p)/n)*exp(seq(from = -8,to = 5,length.out =  100))
  }


  tmp2=VARLASSO(X_t,VAR_order = VAR_order,lambda = lambda,method = "LASSO",crit="bic_m",mbic2_c =mbic2_c)

  if(is.null(bandwidth))
  {
    M0=adaptive_bandwidth_TLM(tmp2$eps,Kern = Kernel,even = T)
    M2=max(4,M0[1])
  }else{
    M2=bandwidth[1]
    M0=c(bandwidth[1],bandwidth[2])
  }

  if(is.null(r))
  {
    N2=M2
  }else{
    N2=floor(M2/(log(M2)^(2/r)))
  }



  FourierFrequencies_h2=1:(N2-1)*pi/N2
  tmp=FourierFrequencies_h2>=Frequencyband[1] & FourierFrequencies_h2<=Frequencyband[2]
  if(sum(tmp)>0)
  {
    FourierFrequencies_h2=FourierFrequencies_h2[tmp]
  }else
  {
    FourierFrequencies_h2=mean(Frequencyband)
  }


  Spec_VAR_esta=spec_inv_preWhiteningVAR_freq(X_t = X_t,lambda=lambda,Kernel = Kernel,M = M2,VAR_order = VAR_order,Frequencies = FourierFrequencies_h2,method = "LASSO",crit="bic_m",Check_pos_def = T,M_corr = M0[2],mbic2_c =mbic2_c,gamma=gamma)

  Test_TX4=max_de_part_coh_H_0_pre_Filt(Z=Z,spec = Spec_VAR_esta$f_inv,Kernel = Kernel,M = M2,DEMEAN = T,Delta_0 = Delta_0,FourierFrequencies_h =  FourierFrequencies_h2,Filter = Array_to_A(Spec_VAR_esta$VAR_A))
  Test_TX_Coh4=max_coh_H_0_pre_Filt(Z=Z,Kernel = Kernel,M = M2,DEMEAN = T,Delta_0 = Delta_0,FourierFrequencies_h =  FourierFrequencies_h2,Filter = Array_to_A(Spec_VAR_esta$VAR_A))

  GM=function(TX,d,Gtx,FDR_alpha)
  {
    t_hat=t_inf(TX[q_ind],Gt=function(t) Gtx(t,d),q = q,ttilde = 2*log(q*(d)),lengt.out = 200,alpha = FDR_alpha)
    return(list(t_hat=t_hat,GM=TX>=t_hat))
  }

  Test1=GM(Test_TX4$TX_chi*n/M2,d=dim(Test_TX4$part_rho)[3],Gtx=Gt2t,FDR_alpha=FDR_alpha)
  Test2=GM(Test_TX_Coh4$TX_chi*n/M2,d=dim(Test_TX_Coh4$coh_rho)[3],Gtx=Gt2t,FDR_alpha=FDR_alpha)



  return(list(Frequencyband_Fourier=FourierFrequencies_h2,Partial_Coherence=Test1$GM,Coherence=Test2$GM))
}

