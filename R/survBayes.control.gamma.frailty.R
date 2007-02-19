`survBayes.control.gamma.frailty` <-
function (mu.cl = 1, prec.tau.cl = 1e-04) 
{
    if (mu.cl < 0) 
        stop("Invalid initialization for mu.cl")
    if (prec.tau.cl < 0) 
        stop("Invalid initialization for prec.tau.cl")
    list(mu.cl = mu.cl, prec.tau.cl = prec.tau.cl, mu.cl.miss = missing(mu.cl), 
        prec.tau.cl.miss = missing(prec.tau.cl))
}

