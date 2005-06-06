\name{nFDR-internal}
\title{Internal nFDR Functions}
\alias{Bf}
\alias{hk}
\alias{f.tilde}
\alias{F.tilde}
\description{
  Internal nFDR functions.
}
\usage{
Bf(f, k, x)
hk(R,k)
F.tilde(tt,x,n)
f.tilde(tt,x,n, Smooth = TRUE)
}

\arguments{
  \item{x}{ a vector of \code{p}-values }
  \item{R}{ vector of \code{r} values }
  \item{f}{ vector of sampled value of function \eqn{f}{f} }
  \item{k}{ order of Bernstein polynomial approximation}
  \item{n}{ order of Bernstein polynomial approximation for cdf \eqn{F}{F} or pdf \eqn{f}{f} }
  \item{tt}{ vector of values on [0,1] at which Bernstein approximation are evaluated }
  \item{Smooth}{ default is TRUE, \eqn{F}{F} is smoothed first }
}

\value{Function \code{bf} returns Bernstein approximation of function \code{f} based on sampled values given in argument \code{f}.
Function \code{hk} returns value of \eqn{h(r,k)}{h(r,k)} used in calculation of the variance of the estimator \eqn{\tilde\pi_0}{tilde_pi_0} of \eqn{\pi_0}{\pi0}. Function \code{F.tilde} returns values of Bernstein approximation \eqn{\tilde F}{F_tilde} of cdf \code{F}. Function \code{f.tilde} returns values of Bernstein approximation \eqn{\tilde f}{f_tilde} of pdf \code{f} with sampled values of \code{f} is numerically estimated by \code{F}, if  \code{Smooth = FALSE}, or by \eqn{\tilde F}{F_tilde}. 
}

\author{ Zhong Guan {zguan@iusb.edu} }

\examples{
}
\keyword{ distribution }
