\name{survBayes}
\alias{survBayes}
\alias{survBayes.b.fctn.Lambda}
\alias{survBayes.base}
\alias{survBayes.chol}
\alias{survBayes.draw.cond.pw.lin.exp}
\alias{survBayes.KM.int}
\alias{survBayes.Lambda0}
\alias{survBayes.numb.events.int}
\alias{survBayes.poisson.update}
\alias{survBayes.taylor}
\title{Fits a proportional hazards model to time to event data by a Bayesian approach}
\description{
Fits a proportional hazards model to time to event data by a Bayesian approach.
Right and interval censored data and a lognormal frailty term can be fitted.
}
\usage{
survBayes(formula = formula(data), data = parent.frame(), burn.in = 1000, number.sample = 1000, max.grid.size = 50, control, control.frailty, seed.set, ...)
}
\arguments{
  \item{formula}{ a formula object, with the response on the left of a \code{~} operator, and the terms on the right.
                 The response must be a survival object of type \code{"right"} or \code{"interval"}
                 as returned by the \code{Surv} function. }
  \item{data}{ a data.frame in which to interpret the variables named in the \code{formula}. }
  \item{burn.in}{ burn.in }
  \item{number.sample}{ number of sample }
  \item{max.grid.size}{ number of grid points }
  \item{control}{ Object of class \code{control} specifying iteration limit and other control options. Default is survBayes.control(...). }
  \item{control.frailty}{ Object of class \code{control.frailty} specifying parameters for the priors of frailties and other control options.
                        Default is survBayes.control.lognormal.frailty(...) or survBayes.control.gamma.frailty(...). }
  \item{seed.set}{ setting of the seed of the random number generator }
  \item{\dots}{ further parameters }
}
\details{
Fits a proportional hazards model to time to event data by a Bayesian approach. The time axis is split into \code{max.grid.size} intervals and the
log baseline hazard is assumed to be a auto regressive process of order one.
Right and interval censored data and a lognormal or gamma frailty term can be fitted. In case of interval censored data the assumed observation times are
augmented by a piecewise exponential distribution conditioned on the respective interval.
}
\value{
The returned values are
  \item{t.where }{used grid points}
  \item{lbh }{samples of the log baseline hazard at the grid points}
  \item{beta}{samples of the vector of covariates}
  \item{sigma.lbh}{samples of sigma.lbh.0 and sigma.lbh.1}
If frailty is included
  \item{alpha.cluster}{samples of the frailty values}
  \item{sigma.cluster}{samples of frailty variance}
The performance of the Metropolis-Hasting steps is checked.
  \item{m.h.performance}{ for beta, lbh and, if appropriate, alpha}
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\seealso{ \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}} }
\examples{
data(AA.data)
AA.res<-survBayes(Surv(t.left,t.right,z*3,type="interval")~mo+lok+frailty(gr,dist="gauss"),data=AA.data,burn.in=10,number.sample=10)
}
\keyword{ survival }
