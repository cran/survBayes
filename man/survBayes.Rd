\name{survBayes}
\alias{survBayes}
\alias{survBayes.b.fctn.Lambda.lin}
\alias{survBayes.b.fctn.Lambda.quad.a}
\alias{survBayes.b.fctn.Lambda.quad.b}
\alias{survBayes.b.fctn.Lambda.quad.c}
\alias{survBayes.b.fctn.Lambda.quad.d}
\alias{survBayes.base}
\alias{survBayes.draw.cond.pw.lin.exp}
\alias{survBayes.int.basis}
\alias{survBayes.KM.int}
\alias{survBayes.poisson.update}
\alias{survBayes.taylor}
\alias{survBayes.time.basis}
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
Fits a proportional hazards model to time to event data by a Bayesian approach. 
The time axis is split into \code{max.grid.size} intervals and the
log baseline hazard is assumed to be cubic spline penalized by an auto regressive process of order one.
Right and interval censored data and a lognormal or gamma frailty term can be fitted. 
In case of interval censored data the assumed observation times are
augmented by a piecewise exponential distribution conditioned on the respective interval.
}
\value{
The returned values are, if appropriate
  \item{t.where }{used grid points}
  \item{beta}{samples of the vector of covariates}
  \item{lbh.coef }{samples of the log baseline hazard coefficients at the grid points}
  \item{sigma.lbh}{samples of sigma.lbh.0 and sigma.lbh.1}
  \item{alpha.cluster}{samples of the frailty values}
  \item{sigma.cluster}{samples of frailty variance}
  \item{z.cluster}{samples of the frailty values}
  \item{mu.cluster}{samples of the rate and shape of the gamma prior}
  \item{m.h.performance}{ The performance of the Metropolis-Hasting steps is checked for beta, lbh and, if appropriate, alpha}
}
\author{ V. Henschel, Ch. Heiss, U. Mansmann }
\seealso{ \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}} }
\examples{
data(AA.data)
control<-survBayes.control(delta.taylor = 0.3, sigma.lbh.1=0.01,rate.sigma.lbh.1 = 1e-3, shape.sigma.lbh.1 = 1e-3)
AA.res<-survBayes(Surv(t.left,t.right,z*3,type="interval")~mo+lok+frailty(gr,dist="gamma"),data=AA.data,burn.in=0,number.sample=10,control=control)
}
\keyword{ survival }
