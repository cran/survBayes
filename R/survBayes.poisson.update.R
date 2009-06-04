"survBayes.poisson.update" <-
function (beta, X.design, cc, offset, R.inv) 
{
    p.par <- dim(X.design)[2]
    a <- rep(0, p.par)
    theta <- X.design %*% beta + offset
    mu <- exp(theta)
    W <- rep(as.vector(mu), p.par)
    y1 <- theta + (cc - mu)/mu - offset
    C <- solve(R.inv + t(X.design) %*% (W * X.design))
    m <- C %*% (R.inv %*% a + t(X.design) %*% (mu * y1))
    LC <- chol(C)
    ind.sample <- rnorm(p.par)
    beta.star <- m + t(LC) %*% ind.sample
    theta.star <- X.design %*% beta.star + offset
    mu.star <- exp(theta.star)
    W.star <- rep(mu.star, p.par)
    y1.star <- theta.star + (cc - mu.star)/mu.star - offset
    C.star <- solve(R.inv + t(X.design) %*% (W.star * X.design))
    m.star <- C.star %*% (R.inv %*% a + t(X.design) %*% (mu.star * 
        y1.star))
    prop <- exp((-t(beta.star - a) %*% R.inv %*% (beta.star - 
        a) + t(beta - a) %*% R.inv %*% (beta - a))/2 + sum(cc * 
        theta.star - cc * theta) + sum(mu - mu.star) + t(beta.star - 
        m) %*% solve(C) %*% (beta.star - m)/2 - t(beta - m.star) %*% 
        solve(C.star) %*% (beta - m.star)/2) * sqrt(det(C)/det(C.star))
    prop <- ifelse(is.na(prop), 0, prop)
    alpha <- min(c(1, prop))
    list(alpha = alpha, beta.star = beta.star, prop = prop)
}
