"survBayes.control.frailty" <-
function (sigma.RE = 100, rate.sigma.clust = 1e-04, shape.sigma.clust = 1e-04) 
{
    if (sigma.RE < 0) 
        stop("Invalid initialization for sigma.RE")
    if (rate.sigma.clust <= 0) 
        stop("Invalid initialization of rate of the prior of rate.sigma.clust")
    if (shape.sigma.clust <= 0) 
        stop("Invalid initialization of shape of the prior of rate.sigma.clust")
    list(sigma.RE = sigma.RE, rate.sigma.clust = rate.sigma.clust, 
        shape.sigma.clust = shape.sigma.clust, sigma.RE.miss = missing(sigma.RE), 
        rate.sigma.clust.miss = missing(rate.sigma.clust), shape.sigma.clust.miss = missing(shape.sigma.clust))
}
