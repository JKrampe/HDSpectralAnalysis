# Transformation Functions ------------------------------------------------
#' X to Z
#'
#'
#' Finite discrete Fourier transform of multivariate time series \eqn{X_t=(X_1,...,X_n)}.
#' Finite discrete Fourier transform \eqn{Z} is to scaling factor \eqn{1/\sqrt(2 \pi n)}.
#'
#' n sample size
#' p dimension
#'
#' @param X_t n times p dimensional matrix. Multivariate time series.
#' @return Z n times p dimensional complex valued matrix. Finite fourer transform.
#' @export
#'
#' @examples
X_to_Z <- function(X_t) {

  n <- dim(X_t)[1]
  return(apply(t(t(X_t) - apply(X_t, 2, mean)) / sqrt(2 * pi * n), 2, fft))
}


#' Z to X
#'
#'
#' Transforms finite discrete Fourier transform \eqn{Z} of multivariate time series \eqn{X_t=(X_1,...,X_n)} back to its multivariate time series
#'
#' n sample size
#' p dimension
#' @param Z n times p dimensional complex valued matrix. Finite fourer transform.
#' @return Z n times p dimensional complex valued matrix. Multivariate time series.
#' @export
#'
#' @examples
Z_to_X=function(Z)
{
  #Finite discrete Fourier transform of X_1,...,X_n
  n=dim(Z)[1]
  return(apply(Z,2,fft,inverse =T)*sqrt(2*pi/n))
}

#' Transform matrix A(p x p*d) to array p times p times d
#'
#' Transform matrix A (p times p*d) to array p times p times d
#' @param A matrix of size p times p*d
#'
#' @return array of size p times p times d
#' @export
#'
#' @examples
A_to_array=function(A)
{
  #Transform matrix A(p times p*d) to array p times p times d
  p=dim(A)[1]
  d=dim(A)[2]/p
  return(sapply(1:d,function(i) A[,1:p+(i-1)*p],simplify = "array"))
}

#' Transform array p times p times d  to matrix A(p times p*d)
#'
#' Transform array p times p times d  to matrix A(p times p*d)
#' @param A array of size p times p times d
#'
#' @return matrix of size p times p*d
#' @export
#'
#' @examples
Array_to_A=function(A)
{
  #Transform array A p times p times d   to matrix A(p times p*d)
  p=dim(A)[1]
  d=dim(A)[3]
  return(matrix(A,nrow = p,ncol = p*d))
}


# Transform Function Spectral ---------------------------------------------

#' f_to_f_inv
#'
#' computes \eqn{f^{-1}} given f for all frequencies of \eqn{f}
#'
#' @param f_0 \eqn{f^{-1}} as array of size p times p times d, where p is the dimension and d is the number of frequencies
#'
#' @return \eqn{f^{-1}} array of size p times p times d, where p is the dimension and d is the number of frequencies
#' @export
#'
#' @examples
f_to_f_inv=function(f_0)
{

  return(array(apply(f_0,3,solve),dim = dim(f_0)))
}

#' f_to_part_rho
#'
#' computes \eqn{\rho} of the partial coherence given \eqn{f} for all frequencies of \eqn{f}.
#'
#' @param f_0 \eqn{f} as array of size p times p times d, where p is the dimension and d is the number of frequencies
#'
#' @return \eqn{\rho}  array of size p times p times d, where p is the dimension and d is the number of frequencies
#' @export
#'
#' @examples
f_to_part_rho=function(f_0)
  {
  #computes rho of partial coherence given f for all frequencies of f
  return(f_inv_to_part_rho(f_to_f_inv(f_0)))
}

f_to_part_coh=function(f_0) #computes partial coherence (|rho|^2) given f for all frequencies of f
  return(abs(f_to_part_rho(f_0))^2)

#' f_inv_to_part_rho
#'
#' computes \eqn{\rho} of the partial coherence given \eqn{f^{-1}} for all frequencies of \eqn{f^{-1}}.
#'
#' @param f_0 \eqn{f^{-1}} as array of size p times p times d, where p is the dimension and d is the number of frequencies
#'
#' @return \eqn{\rho}  array of size p times p times d, where p is the dimension and d is the number of frequencies
#' @export
#'
#' @examples
f_inv_to_part_rho=function(f_0)
{
  #
  p=dim(f_0)[1]
  return(array(apply(f_0,3,function(x) diag(1/sqrt(Re(diag(x))))%*%(-x)%*%diag(1/sqrt(Re(diag(x))))+diag(p)*2),dim = dim(f_0)))
}

#' f_to_Gamma
#'
#' Computes autocovariance \eqn{\Gamma} to lag h given spectral denstiy f_0
#'
#' @param f_0 \eqn{f} as array of size p times p times d, where p is the dimension and d is the number of frequencies
#' @param h lag for the autocovariance
#'
#' @return \eqn{\Gamma} as array of size p times p
#' @export
#'
#' @examples
f_to_Gamma=function(f_0,h)
{
  #Computes Autocovariance to lag h given spectral denstiy f_0
  Fouriern=dim(f_0)[3]
  FourierFrequencies=(0:(Fouriern-1))/Fouriern*(2*pi)
  return(apply(aperm(f_0,c(3,1,2))*exp(1i*FourierFrequencies*h),2:3,mean))
}

#' Coherence
#'
#' computes the coherence given \eqn{f} for all frequencies of \eqn{f}.
#'
#' @param f \eqn{f} as array of size p times p times d, where p is the dimension and d is the number of frequencies
#'
#' @return coherence as array of size p times p times d, where p is the dimension and d is the number of frequencies
#' @export
#'
#' @examples
f_to_norm_f=function(f)
{
  #Computes s of coherence given f for all frequences of f
  f_diag=1/sqrt(Re(apply(f,3,diag)))
  return(sapply(1:dim(f)[3],function(i) diag(f_diag[,i])%*%f[,,i]%*%diag(f_diag[,i]),simplify = "array"))
}

#' Gamma_to_f
#'
#' Computes the spectral density matrix at frequency omega given an array of autocovariance matrices (lag=0,1,...,d-1).
#'
#' @param Gamma_dn_array array of size p times p times d, where p is the dimension and d is the number of autocovariance. The lags need to be stored as lag=0,1,...,d-1
#' @param omega frequency, where the spectral density matrix should be computed
#'
#' @return spectral density matrix at frequency omega as p times p matrix.
#' @export
#'
#' @examples
Gamma_to_f=function(Gamma_dn_array,omega)
{
  #Given array of ACF computes spectral density at frequencies omega
  p=dim(Gamma_dn_array)[3]
  Gamma_dn_array2=sapply(p:2,function(i) t(Gamma_dn_array[,,i]),simplify = "array")
  Gamma_dn_array2=abind::abind(Gamma_dn_array2,Gamma_dn_array,along = 3)
  return(apply(sapply((-p+1):(p-1),function(i) Gamma_dn_array2[,,i+p]*exp(-1i*(i)*omega),simplify = "array"),1:2,sum))
}
