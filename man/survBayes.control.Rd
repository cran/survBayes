\name{survBayes.control}
\alias{survBayes.control}
\title{ Package options for survBayes }
\description{
These function checks and packages the fitting options for \code{survBayes}
}
\usage{
survBayes.control(n.inter = 10, delta.taylor = 0.1, haz.global = Inf, sigma.lbh.0 = 100, sigma.lbh.1 = 100, 
prec.beta = 1e-04, rate.sigma.lbh.0 = 1e-04, rate.sigma.lbh.1 = 1e-04, shape.sigma.lbh.0 = 1e-04, shape.sigma.lbh.1 = 1e-04,
beta.init=NULL)
}
\arguments{
  \item{n.inter}{ number of intervals to display}
  \item{delta.taylor}{ bandwidth for Taylor approximation }
  \item{haz.global}{ initialization of the global hazard }
  \item{sigma.lbh.0}{ initialization of \code{sigma.lbh.0} }
  \item{sigma.lbh.1}{ initialization of \code{sigma.lbh.0} }
  \item{prec.beta}{ initialization of the precision of the prior of \code{beta}}
  \item{rate.sigma.lbh.0}{ initialization of rate of the prior of \code{sigma.lbh.0} }
  \item{rate.sigma.lbh.1}{ initialization of rate of the prior of \code{sigma.lbh.1} }
  \item{shape.sigma.lbh.0}{ initialization of shape of the prior of \code{sigma.lbh.0} }
  \item{shape.sigma.lbh.1}{ initialization of shape of the prior of \code{sigma.lbh.1} }
  \item{beta.init}{ initialization of \code{beta} }
}
\value{
  A list with the same elements as the input
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\examples{

}
\keyword{ survival }
