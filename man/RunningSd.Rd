% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convo_func.R
\name{RunningSd}
\alias{RunningSd}
\title{Fast Running Standard Deviation Computation}
\usage{
RunningSd(x, W, circular = FALSE)
}
\arguments{
\item{x}{A numeric vector.}

\item{W}{A numeric scalar; width of \code{x} window over which sample variance is computed.}

\item{circular}{Logical; whether  running sample standard deviation is computed assuming
circular nature of  \code{x} sequence (see Details).}
}
\value{
A numeric vector.
}
\description{
Computes running sample standard deviation of a sequence in a fixed width window. Uses
convolution implementation via Fast Fourier Transform.
}
\details{
The length of output vector equals the length of \code{x} vector.
Parameter \code{circular} determines whether \code{x} sequence is assumed to have a  circular nature.
Assume \eqn{l_x} is the length of sequence \code{x}, \code{W} is a fixed length of \code{x} sequence window.

If \code{circular} equals \code{TRUE} then
  \itemize{
    \item first element of the output sequence corresponds to sample standard deviation of \code{x[1:W]},
    \item last element of the output sequence corresponds to sample standard deviation of \code{c(x[l_x], x[1:(W - 1)])}.
  }

If \code{circular} equals \code{FALSE} then
  \itemize{
    \item first element of the output sequence corresponds to sample standard deviation of \code{x[1:W]},
    \item the \eqn{l_x - W + 1}-th element of the output sequence corresponds to sample standard deviation of \code{x[(l_x - W + 1):l_x]},
    \item last \code{W-1} elements of the output sequence are filled with \code{NA}.
  }
}
\examples{
x <- rnorm(10)
RunningSd(x, 3, circular = FALSE)
RunningSd(x, 3, circular = FALSE)

}
