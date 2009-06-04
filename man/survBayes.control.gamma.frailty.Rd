\name{survBayes.control.gamma.frailty}
\alias{survBayes.control.gamma.frailty}
\title{ Package options for gamma frailty in survBayes }
\description{
  These function checks and packages the fitting options for the gamma frailty in \code{survBayes}.
  The frailty values are assumed the be gamma distributed with rate and shape \code{mu.cl} such
  that the expected value is one. The prior of tau.cl=log(mu.cl) is assumed to be normal distributed 
  with mean zero and precision \code{prec.tau.cl}
}
\usage{
survBayes.control.gamma.frailty(mu.cl = 1, prec.tau.cl = 1e-04)
}
\arguments{
  \item{mu.cl}{ initialization of rate and shape of the prior of \code{mu.cl} }
  \item{prec.tau.cl}{ initialization of precision of the prior of \code{prec.tau.cl} }
}
\value{
  A list with the same element as the input
}
\author{ V. Henschel, U. Mansmann }
\keyword{ survival }
