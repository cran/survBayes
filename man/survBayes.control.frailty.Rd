\name{survBayes.control.frailty}
\alias{survBayes.control.frailty}
\title{ Package options for frailty in survBayes }
\description{
  These function checks and packages the fitting options for the gaussian frailty in \code{survBayes}
}
\usage{
survBayes.control.frailty(sigma.RE = 100, rate.sigma.clust = 1e-04, shape.sigma.clust = 1e-04)
}
\arguments{
  \item{sigma.RE}{ initialization of \code{sigma.RE} }
  \item{rate.sigma.clust}{ initialization of rate of the prior of \code{rate.sigma.clust} }
  \item{shape.sigma.clust}{ initialization of shape of the prior of \code{shape.sigma.clust} }
}
\value{
  A list with the same elements as the input
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\examples{

}
\keyword{ survival }
