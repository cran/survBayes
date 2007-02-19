\name{survBayes.baseline.hazard}
\alias{survBayes.baseline.hazard}
\title{ Baseline hazard of survBayes result }
\description{
These function calculates the plain, log or cumulative baseline hazard for a \code{survBayes} result
}
\usage{
survBayes.baseline.hazard(surv.res, type = "log", ci = FALSE, n.inter = 3, start = NULL, end = NULL, thin = 1)
}
\arguments{
  \item{surv.res}{ result of \code{survBayes} }
  \item{type}{ One of \code{"log"}(default),\code{"plain"},\code{"cum"}. Determines if the 
        log baseline hazard, the baseline hazard or the cumulative baseline hazard is calculated.}
  \item{ci}{ if TRUE credibility intervals of the chosen type of baseline hazard are calculated}
  \item{n.inter}{ number of points between the interval points to display, not for type="cum" }
  \item{start}{ the first iteration of interest }
  \item{end}{ the last iteration of interest }
  \item{thin}{ the required interval between successive samples }
}
\value{
The returned values are
  \item{ time }{used times}
  \item{log.base.haz}{ log baseline hazard, if type=\code{"log"}, with \code{.lower} and \code{.upper} if ci }
  \item{base.haz}{ baseline hazard, if type=\code{"plain"}, with \code{.lower} and \code{.upper} if ci }
  \item{cum.base.haz}{ cumulative baseline hazard, if type=\code{"cum"}, with \code{.lower} and \code{.upper} if ci }
}
\author{ V. Henschel, U. Mansmann }
\seealso{ \code{\link[coda]{window.mcmc}} }
\examples{

}
\keyword{ survival }
