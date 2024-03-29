% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spectral.R
\name{var_part_coherence_inv}
\alias{var_part_coherence_inv}
\title{Chi^2 version of de-biased partial coherence}
\usage{
var_part_coherence_inv(x, y)
}
\arguments{
\item{x}{vectorized array of partial coherences of size p times p times d}

\item{y}{vectorized array of test statistics of size p times p times d}
}
\value{
Chi^2 version of de-biased partial coherence
}
\description{
Computes the Chi^2 version of de-biased partial coherence. For this, first the variance matrix of the real and imaginary part is computed.
Then, the inverse matrix of this is used as a transformtion on the real and imaginary parts.
}
