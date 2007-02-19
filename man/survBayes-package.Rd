\name{survBayes-package}
\alias{survBayes-package}
\docType{package}
\title{
Fits a proportional hazards model to time to event data by a Bayesian approach}
\description{
Fits a proportional hazards model to time to event data by a Bayesian approach.
Right and interval censored data and a lognormal or gamma frailty term can be fitted.}
\details{
\tabular{ll}{
Package: \tab survBayes\cr
Type: \tab Package\cr
Version: \tab 0.2.1\cr
Date: \tab 2007-02-19\cr
License: \tab GPL Version 2 or newer\cr
}
Fits a proportional hazards model to time to event data by a Bayesian approach. 
The time axis is split into \code{max.grid.size} intervals and the
log baseline hazard is assumed to be cubic spline penalized by an auto regressive process of order one.
Right and interval censored data and a lognormal or gamma frailty term can be fitted. 
In case of interval censored data the assumed observation times are
augmented by a piecewise exponential distribution conditioned on the respective interval.
}
\author{
Volkmar Henschel, Christiane Heiss, Ulrich Mansmann

Maintainer: Volkmar Henschel <henschel@ibe.med.uni-muenchen.de>
}
\keyword{ package }
\keyword{ survival }
\seealso{ \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}} }
\examples{
data(aneurism.data)
control<-survBayes.control(delta.taylor = 0.3, sigma.lbh.1=0.01,rate.sigma.lbh.1 = 1e-3, shape.sigma.lbh.1 = 1e-3)
aneurism.res<-survBayes(Surv(left,right,cens*3,type="interval")~mo+loc+frailty(gr,dist="gamma"),data=aneurism.data,burn.in=0,number.sample=10,control=control)
}
