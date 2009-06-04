"survBayes.control" <-
function (n.inter = 100, delta.taylor = 0.1, sigma.lbh.0 = 100, 
    sigma.lbh.1 = 100, prec.beta.init = 1e-04, rate.wishart.beta = 1e-04, 
    shape.wishart.beta = 1e-04, rate.sigma.lbh.0 = 1e-04, rate.sigma.lbh.1 = 1e-04, 
    shape.sigma.lbh.0 = 1e-04, shape.sigma.lbh.1 = 1e-04, beta.init = NULL) 
{
    if (n.inter < 0) 
        stop("Invalid value for n.inter")
    if (delta.taylor <= 0) 
        stop("Invalid width for L2 approximation of exp")
    if (sigma.lbh.0 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (sigma.lbh.1 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (prec.beta.init <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (rate.sigma.lbh.0 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (rate.wishart.beta <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (shape.wishart.beta <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (rate.sigma.lbh.1 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (shape.sigma.lbh.0 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    if (shape.sigma.lbh.1 <= 0) 
        stop("Invalid parameters/values for prior distributions")
    list(n.inter = n.inter, delta.taylor = delta.taylor, sigma.lbh.0 = sigma.lbh.0, 
        sigma.lbh.1 = sigma.lbh.1, prec.beta.init = prec.beta.init, 
        rate.wishart.beta = rate.wishart.beta, shape.wishart.beta = shape.wishart.beta, 
        rate.sigma.lbh.0 = rate.sigma.lbh.0, rate.sigma.lbh.1 = rate.sigma.lbh.1, 
        shape.sigma.lbh.0 = shape.sigma.lbh.0, shape.sigma.lbh.1 = shape.sigma.lbh.1, 
        beta.init = beta.init, n.inter.miss = missing(n.inter), 
        delta.taylor.miss = missing(delta.taylor), sigma.lbh.0.miss = missing(sigma.lbh.0), 
        sigma.lbh.1.miss = missing(sigma.lbh.1), prec.beta.init.miss = missing(prec.beta.init), 
        rate.wishart.beta.miss = missing(rate.wishart.beta), 
        shape.wishart.beta.miss = missing(shape.wishart.beta), 
        rate.sigma.lbh.0.miss = missing(rate.sigma.lbh.0), rate.sigma.lbh.1.miss = missing(rate.sigma.lbh.1), 
        shape.sigma.lbh.0.miss = missing(shape.sigma.lbh.0), 
        shape.sigma.lbh.1.miss = missing(shape.sigma.lbh.1))
}
