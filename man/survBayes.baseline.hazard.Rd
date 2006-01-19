\name{survBayes.baseline.hazard}
\alias{survBayes.baseline.hazard}
\title{ Baseline hazard of survBayes result }
\description{
These function calculates the plain, log or cumulative baseline hazard for a \code{survBayes} result
}
\usage{
survBayes.baseline.hazard(surv.res,type="log",n.inter=3,start=NULL,end=NULL,thin=1)
}
\arguments{
  \item{surv.res}{ result of \code{survBayes} }
  \item{type}{ One of \code{"log"}(default),\code{"plain"},\code{"cum"}. Determines if the 
        log baseline hazard, the baseline hazard or the cumulative baseline hazard is calculated.}
  \item{n.inter}{ number of points between the interval points to display }
  \item{start}{ the first iteration of interest }
  \item{end}{ the last iteration of interest }
  \item{thin}{ the required interval between successive samples }
}
\value{
The returned values are
  \item{ time }{used times}
  \item{log.base.haz}{ log baseline hazard, if type=\code{"log"} }
  \item{base.haz}{ baseline hazard, if type=\code{"plain"} }
  \item{cum.base.haz}{ cumulative baseline hazard, if type=\code{"cum"} }
}
\author{ V. Henschel, U. Mansmann }
\seealso{ \code{\link[coda]{window.mcmc}} }
\examples{

}
\keyword{ survival }
