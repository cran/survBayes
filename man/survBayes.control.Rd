\name{survBayes.control}
\alias{survBayes.control}
\title{ Package options for survBayes }
\description{
These function checks and packages the fitting options for \code{survBayes}
}
\usage{
survBayes.control(n.inter = 1000, delta.taylor = 0.3, beta.init = NULL, prec.beta.init = 1e-04, 
sigma.lbh.0 = 100, sigma.lbh.1 = 100, 
rate.sigma.lbh.0 = 1e-04, rate.sigma.lbh.1 = 1e-04, 
shape.sigma.lbh.0 = 1e-04, shape.sigma.lbh.1 = 1e-04)
}
\arguments{
  \item{n.inter}{ number of intervals to display}
  \item{delta.taylor}{ bandwidth for Taylor approximation }
  \item{beta.init}{ initialization of \code{beta} }
  \item{prec.beta.init}{ initialization of the precision of the prior of \code{beta} }
  \item{sigma.lbh.0}{ initialization of \code{sigma.lbh.0} }
  \item{sigma.lbh.1}{ initialization of \code{sigma.lbh.0} }
  \item{rate.sigma.lbh.0}{ initialization of rate of the gamma prior of \code{sigma.lbh.0} }
  \item{rate.sigma.lbh.1}{ initialization of rate of the gamma prior of \code{sigma.lbh.1} }
  \item{shape.sigma.lbh.0}{ initialization of shape of the gamma prior of \code{sigma.lbh.0} }
  \item{shape.sigma.lbh.1}{ initialization of shape of the gamma prior of \code{sigma.lbh.1} }
}
\value{
  A list with the same elements as the input
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\examples{

}
\keyword{ survival }
