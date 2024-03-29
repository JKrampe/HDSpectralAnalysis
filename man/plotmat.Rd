% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{plotmat}
\alias{plotmat}
\title{Plot a matrix}
\usage{
plotmat(
  A,
  color1 = "red",
  color2 = "blue",
  range_A = NULL,
  xlab = "column",
  ylab = "row",
  pch = 15,
  alpha_func = function(x, y) scales::alpha(x, y),
  ...
)
}
\arguments{
\item{A}{matrix to plot}

\item{color1}{Color for positive entries}

\item{color2}{Color for negative entries}

\item{range_A}{Value range to consider}

\item{xlab}{Name for x-axis}

\item{ylab}{name for y-axis}

\item{pch}{pch-number used for plotted points (15 default)}

\item{alpha_func}{Function used for color-transparency}

\item{...}{Additional parameters for the plot function}
}
\value{
No return value
}
\description{
This functions plots a matrix A, i.e, x- axis are the column indices, y-axis are the row indices.
Color1 is used for positive, Color2 for negative entries. Zero entries are left blank. Entries closer to zero are plotted more transparent.
}
