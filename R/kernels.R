# Kernels -----------------------------------------------------------------
#with compact support in time domain written in frequency domain



#' Uniform kernel
#'
#' Uniform kernel with lag window M written in frequency domain.
#' Note, this kernel is not positive definite.
#'
#' @param lambda input frequencies
#' @param M window length
#'
#' @return kernel weights
#' @export
#'
#' @examples
Uniform=function(lambda,M)
{
  return((sin((M+1/2)*lambda)/(.Machine$double.eps+sin(lambda/2)))+as.numeric(lambda%in%c(0,-2*pi,2*pi))*(2*M+1))
}



#' Modified Bartlett kernel
#'
#' Modified Bartlett kernel (triangular kernel \eqn{K(u)=1-|u|} for \eqn{u \in (-1,1)}) with lag window M written in frequency domain.
#' Note, this kernel is positive definite.
#'
#' @param lambda input frequencies
#' @param M window length
#'
#' @return kernel weights
#' @export
#'
#' @examples
Modified_Bartlett=function(lambda,M)
{
  return(1/(M)*(sin(M*lambda/2)^2/(.Machine$double.eps+sin(lambda/2)^2)+as.numeric(lambda%in%c(0,-2*pi,2*pi))*(M^2)))
}
