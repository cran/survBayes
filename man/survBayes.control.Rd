\name{survBayes.control}
\alias{survBayes.control}
\title{ Package options for survBayes }
\description{
These function checks and packages the fitting options for \code{survBayes}
}
\usage{
survBayes.control(n.inter = 100, delta.taylor = 0.1, sigma.lbh.0 = 100, sigma.lbh.1 = 100, 
prec.beta.init = 1e-04, rate.wishart.beta = 1e-04, shape.wishart.beta = 1e-04, 
rate.sigma.lbh.0 = 1e-04, rate.sigma.lbh.1 = 1e-04, shape.sigma.lbh.0 = 1e-04, shape.sigma.lbh.1 = 1e-04,
beta.init=NULL)
}
\arguments{
  \item{n.inter}{ number of intervals to display}
  \item{delta.taylor}{ bandwidth for Taylor approximation }
  \item{sigma.lbh.0}{ initialization of \code{sigma.lbh.0} }
  \item{sigma.lbh.1}{ initialization of \code{sigma.lbh.0} }
  \item{prec.beta.init}{ initialization of the precision of the prior of \code{beta} }
  \item{rate.wishart.beta}{ initialization of rate of the wishart prior of \code{cov.beta} }
  \item{shape.wishart.beta}{ initialization of shape factor of the wishart prior of \code{cov.beta} 
        the matrix is 1 at the diagonal and 0.05 at the off diagonal elements}
  \item{rate.sigma.lbh.0}{ initialization of rate of the gamma prior of \code{sigma.lbh.0} }
  \item{rate.sigma.lbh.1}{ initialization of rate of the gamma prior of \code{sigma.lbh.1} }
  \item{shape.sigma.lbh.0}{ initialization of shape of the gamma prior of \code{sigma.lbh.0} }
  \item{shape.sigma.lbh.1}{ initialization of shape of the gamma prior of \code{sigma.lbh.1} }
  \item{beta.init}{ initialization of \code{beta} }
}
\value{
  A list with the same elements as the input
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\keyword{ survival }
