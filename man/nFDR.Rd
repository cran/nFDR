\name{nFDR}
\alias{nFDR}
\title{ Nonparametric Estimate of FDR Based on Bernstein Polynomials }
\description{ Using the Bernstein polynomial
density estimation to estimate the proportion \eqn{\pi_0}{\pi0} of true null hypotheses based \eqn{p}{p}-values. 
 A nonparametric estimate of false discovery
rates (FDR's) are also calculated.
  \code{nFDR} calls a C program to search the minimizer \eqn{(r,k)}{(r,k)} of a partial mean square error and returns the proportion of true null hypotheses, FDR, FNR, and q-values.
}

\usage{
nFDR(x, r0, k0, K, alpha = 0.05, Trial.r = 5, Trial.k = 10, Method = "approx", Smooth = TRUE)
}

\arguments{
  \item{x}{ a vector of p-values }
  \item{r0}{ initial guess of r }
  \item{k0}{ initial guess of k }
  \item{K}{ Upper limit for r }
  \item{alpha}{ significance level }
  \item{Trial.r}{ number of trials for searching r }
  \item{Trial.k}{ number of trials for searching k }
  \item{Method}{ default is "approx" }
  \item{Smooth}{ default is TRUE }
}
\details{If either of \code{Trial.r} or \code{Trial.k} is zero, then the search of \eqn{(r,k)}{(r,k)} is skipped and choose \eqn{(r,k)=(r_0,k_0)}{(r,k)=(r0,k0)}.
}

\value{
  a list containing: 
  \item{p}{ sorted p-values }
  \item{(r,k) }{minimizer of the partial mean square error pMSE }
  \item{PI0 }{ estimated proportion \eqn{\pi_0}{\pi0} of true null hypotheses }
  \item{cint.pi0 }{ confidence interval for the proportion \eqn{\pi_0}{\pi0} of true null hypotheses }
  \item{FDR }{ estimated false discovery rate (FDR) according to \eqn{p}{p}}
  \item{FNR }{ estimated false nondiscovery rate (FNR) according to \eqn{p}{p}}
  \item{Cint.fdr }{ confidence interval for FDR according to \eqn{p}{p}}
  \item{qvalue }{ qvalue according to \eqn{p}{p}}
}

\author{ Zhong Guan {zguan@iusb.edu} }
\references{
Zhong Guan, Baolin Wu and Hongyu Zhao(2005),   
	Nonparametric estimator of false discovery rate based on Bernstein polynomials  
}

\examples{
set.seed (777)
library(nFDR) # load the package
n<-200  # sample size
pi0<-0.70 # proportion of true null hypotheses
x<-c(runif(n*pi0,0,1), rbeta(n*(1-pi0), 1, 6))
 ## simulate n p-values from mixture of beta(1,6) and uniform(0,1)
res<-nFDR(x, r0 = 5, k0 = 50, K = .5*n, alpha = 0.05, Trial.r = 3, Trial.k = 5, Method = "approx", Smooth = TRUE)
res$PI0 # estimated proportion of true nulls
}
\keyword{ htest}
\keyword{ nonparametric}
\keyword{ smooth }
