`survBayes.base` <-
function (int.matrix, type.ind, X.design, frailty.values, frailty.dist, 
    formula, burn.in, number.sample, max.grid.size, data, control, 
    control.frailty, seed.set) 
{
    set.seed(seed.set)
    n.inter <- control$n.inter                                          # Interval to display number of cycles
    delta.taylor <- control$delta.taylor                                # Width for L2 approximation of exp
#
#       Setting up the data
#
#       In case of right censored data, create time and censor variable in the classical form
#
    if (type.ind == 1) {
        time <- int.matrix[, 1]
        time.max <- max(int.matrix[, 1])
        cens <- int.matrix[, 2]
    }
#
#       In case of finite interval censored data create time as the midpoint of the interval
#       Keep the right censred observations in the classical form 
#
    else {
        time <- ifelse(int.matrix[, 3] == 3, (int.matrix[, 1] + 
            int.matrix[, 2])/2, int.matrix[, 1])
        time.max <- max(c(int.matrix[, 1], int.matrix[, 2]), 
            na.rm = TRUE)
        cens <- int.matrix[, 3]/3
    }
#
#
    if (!is.matrix(X.design)) 
        X.design <- matrix(X.design, ncol = 1)
    pp <- ncol(X.design)                                                # number of covariates
    nn <- nrow(X.design)                                                # number of observations
#
#       Setting up the time grid, time intervals
#       int.left are the left boundries of the intervals
#       int.right are the right boundries of the intervals
#       int.delta is the length of the single intervals
#       max.grid.size is the number of intervals
#       int are the "quasi quantiles" calculated from a linearized KM curve
#       survBayes.KM.int(time, cens, NN) creates NN-1 intervals determined by NN time points
#
    KK <- max.grid.size + 2
    int <- survBayes.KM.int(time, cens, max.grid.size + 1)
    int[length(int)] <- time.max
    int.left <- int[-length(int)]
    int.right <- int[-1]
    int.delta <- diff(int)
#
#       The use of spline motivates to add to additional time points before the first int element
#       as well as two time points after the last int element. 
#
    int.help <- c(int[1] - 2 * int.delta[1], int[1] - int.delta[1], 
        int, time.max + 1, time.max + 2)
    int.delta.help <- apply((cbind(diff(int.help)[-c(KK + 1, 
        KK + 2)], diff(int.help)[-c(1, KK + 2)], diff(int.help)[-c(1, 
        2)])), 1, mean)
    BB.int <- sapply(1:(KK - 1), survBayes.int.basis, int.help,         #       Matrix of the Basis functions at the interval points
        KK)
#
#       Inits needed:
#       Inits coefficients (beta)
#       Inits log baseline hazard (lbh)
#       lbh is initialized by taking for h the mean incidence per time interval,
#               the last interval gets lbh=0
#       Inits variability for lbh 
#                       sigma.lbh.0 - variance for starting point
#                       sigma.lbh.1 - variance for smoothness of walk
#
    beta.init <- control$beta.init
    if (is.null(beta.init)) {
        beta <- rep(0, pp)
    }
    else {
        if (length(beta.init) == pp) {
            beta <- beta.init
        }
        else {
            stop("Invalid length of beta.init")
        }
    }
    index.vec <- 1:(KK + 1)
    lbh.coef <- rep(0, KK + 1)
#
#       Parameters for prior distributions needed
#
#        haz.global <- sum(cens)/sum(time)
#        lbh.coef <- rep(log(haz.global), KK +1)
#
    sigma.lbh.0 <- control$sigma.lbh.0
    sigma.lbh.1 <- control$sigma.lbh.1
    prec.beta.init <- control$prec.beta.init
    Cov.beta.inv <- diag(pp) * prec.beta.init
    rate.sigma.lbh.0 <- control$rate.sigma.lbh.0
    rate.sigma.lbh.1 <- control$rate.sigma.lbh.1
    shape.sigma.lbh.0 <- control$shape.sigma.lbh.0
    shape.sigma.lbh.1 <- control$shape.sigma.lbh.1
#
#       Precision matrix components for random walk
#       The following sets the framework for a first order random walk (QQ.ar1)
#       as well as for a second order random walk (QQ.ar2).
#       For second order random walk we also need a dispersion parameter for the second point 
#       of the random walk.
#
    QQ.0 <- matrix(0, nrow = KK + 1, ncol = KK + 1)
    QQ.0[1, 1] <- 1
    QQ.1 <- matrix(0, nrow = KK + 1, ncol = KK + 1)
    QQ.1[2, 2] <- 1
    QQ.ar1 <- diag(c(1/int.delta.help, 0) + c(0, 1/int.delta.help))
    for (i in 2:(KK + 1)) {
        QQ.ar1[i - 1, i] <- (-1)/int.delta.help[i - 1]
        QQ.ar1[i, i - 1] <- (-1)/int.delta.help[i - 1]
    }
    QQ.ar2 <- diag(c(1, 5, rep(6, KK - 3), 5, 1))
    for (i in 3:KK) {
        QQ.ar2[i - 1, i] <- -4
        QQ.ar2[i, i - 1] <- -4
    }
    QQ.ar2[1, 2] <- -2
    QQ.ar2[2, 1] <- -2
    QQ.ar2[KK + 1, KK] <- -2
    QQ.ar2[KK, KK + 1] <- -2
    for (i in 3:(KK + 1)) {
        QQ.ar2[i - 2, i] <- 1
        QQ.ar2[i, i - 2] <- 1
    }
#
#       the frailty term is set
#
    if (is.null(frailty.values)) {
        frailty.coef <- rep(0, nn)
    }
    else {
        n.cl <- length(unique(frailty.values))
        if (frailty.dist == "gauss") {
#       Inits needed:
#       alpha.cluster contains the cluster specific random effect estimates
#       sigma.RE is the variance of the randomeffect distribution N(0,sigma.re)
            alpha.cluster <- rep(0, n.cl)
            sigma.RE <- control.frailty$sigma.RE
            rate.sigma.clust <- control.frailty$rate.sigma.clust
            shape.sigma.clust <- control.frailty$shape.sigma.clust
            alpha.cluster.out <- NULL
            sigma.RE.out <- NULL
            frailty.coef <- alpha.cluster[frailty.values]
        }
        else {
            if (frailty.dist == "gamma") {
#           Inits needed:
#           z.cluster contains the cluster specific random effect estimates
#           mu.cl is the rate and the shape of the random effect distribution Ga(mu.cl,mu.cl)
                z.cluster <- rep(1, n.cl)
                mu.cl <- control.frailty$mu.cl
                tau.cl <- log(mu.cl)
                prec.tau.cl <- 1e-04
                z.cluster.out <- NULL
                mu.cl.out <- NULL
                frailty.coef <- log(z.cluster[frailty.values])
            }
            else {
                stop("Invalid frailty distribution")
            }
        }
    }
#
#       Preparing the sampling
#
    iter <- 0
    beta.out <- NULL
    lbh.coef.out <- NULL
    sigma.lbh.out <- NULL
    m.h.performance <- NULL
#
#       Starting the sampling
#
    sample.total <- burn.in + number.sample
    while (iter < sample.total) {
        iter <- iter + 1
        m.h.beta <- 0
        m.h.lbh <- 0
        m.h.alpha <- 0
        m.h.cl <- 0
#
#       make the data augmentation
#
        lbh.int <- as.vector(lbh.coef %*% BB.int)
        Lambda0.int <- int.delta/2 * (exp(lbh.int[-(KK - 1)]) + 
            exp(lbh.int[-1]))
        Lambda0 <- c(0, cumsum(Lambda0.int))
        log.hr <- frailty.coef + X.design %*% beta
        if (type.ind == 2) {
            inf.pat.mat <- cbind(int.matrix[, 1:2], log.hr)
            time.aug <- apply(inf.pat.mat, 1, survBayes.draw.cond.pw.lin.exp, 
                int = int, Lambda0 = Lambda0)
            time <- time.aug[1, ]
            cens <- time.aug[2, ]
        }
#
#               Perform the beta update
#
        time.int <- sapply(time, function(t, int.left, int.right) {
            which(int <= t & c(int.right, time.max + 1) > t)
        }, int.left, int.right)
        time.int.mat <- t(sapply(time, function(t, int.right) {
            as.numeric(int.right <= t)
        }, int.right))
        BB.time <- sapply(time, survBayes.time.basis, int.help, 
            KK)
        Lambda0.int.last <- (time - int[time.int])/2 * (exp(lbh.int[time.int]) + 
            exp(as.vector(lbh.coef %*% BB.time)))
        Lambda0 <- Lambda0[time.int] + Lambda0.int.last
        offset <- frailty.coef + log(Lambda0)
        coef.update <- survBayes.poisson.update(beta, X.design, 
            cens, offset, Cov.beta.inv)
        u <- runif(1)
        if (u <= coef.update$alpha) {
            beta <- coef.update$beta.star
            m.h.beta <- 1
        }
        if (iter > burn.in) 
            beta.out <- cbind(beta.out, beta)
#
#               Perform the lbh update following section 3.1.2 of Rue (2001)
#
#               The b-term in fomula (4) of Rue 2001 is calculated
#
        hr <- exp(frailty.coef + X.design %*% beta)
        taylor.coef.int <- survBayes.taylor(as.vector(lbh.coef %*% 
            BB.int), delta.taylor)
        taylor.coef.time <- survBayes.taylor(as.vector(lbh.coef %*% 
            BB.time), delta.taylor)
        b.vec.cens <- as.vector(cens %*% t(BB.time))
        b.vec.Lambda.lin <- sapply(index.vec, survBayes.b.fctn.Lambda.lin, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int$b, 
            taylor.time = taylor.coef.time$b)
        b.vec.Lambda.quad.a <- diag(sapply(index.vec, survBayes.b.fctn.Lambda.quad.a, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int$c, 
            taylor.time = taylor.coef.time$c))
        b.vec.Lambda.quad.b <- cbind(rep(0, KK + 1), diag(sapply(1:KK, 
            survBayes.b.fctn.Lambda.quad.b, time = time, time.int = time.int, 
            time.int.mat = time.int.mat, hr = hr, int = int, 
            int.delta = int.delta, BB.int = BB.int, BB.time = BB.time, 
            KK = KK, taylor.int = taylor.coef.int$c, taylor.time = taylor.coef.time$c), 
            nrow = KK + 1, ncol = KK))
        b.vec.Lambda.quad.c <- cbind(rep(0, KK + 1), rep(0, KK + 
            1), diag(sapply(1:(KK - 1), survBayes.b.fctn.Lambda.quad.c, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int$c, 
            taylor.time = taylor.coef.time$c), nrow = KK + 1, 
            ncol = KK - 1))
        b.vec.Lambda.quad.d <- cbind(rep(0, KK + 1), rep(0, KK + 
            1), rep(0, KK + 1), diag(sapply(1:(KK - 2), survBayes.b.fctn.Lambda.quad.d, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int$c, 
            taylor.time = taylor.coef.time$c), nrow = KK + 1, 
            ncol = KK - 2))
        b.vec.Lambda.quad <- b.vec.Lambda.quad.a + b.vec.Lambda.quad.b + 
            t(b.vec.Lambda.quad.b) + b.vec.Lambda.quad.c + t(b.vec.Lambda.quad.c) + 
            b.vec.Lambda.quad.d + t(b.vec.Lambda.quad.d)
        b.vec <- b.vec.cens - b.vec.Lambda.lin
#
#               The Q-term in fomula (4) of Rue (2001) is calculated
#
        QQ.hlp <- QQ.0/sigma.lbh.0 + QQ.ar1/sigma.lbh.1
        QQ <- QQ.hlp + b.vec.Lambda.quad
#
#               Sampling a new candidate for lbh following section 3.1.2 of Rue (2001)
#
        LL <- t(chol(QQ))
        mu.lbh <- solve(t(LL), solve(LL, b.vec))
        lbh.coef.star <- mu.lbh + solve(t(LL), rnorm(KK + 1))
#
#               Metropolis Hastings step to accept the new candidate
#
        taylor.coef.int.star <- survBayes.taylor(as.vector(lbh.coef.star %*% 
            BB.int), delta.taylor)
        taylor.coef.time.star <- survBayes.taylor(as.vector(lbh.coef.star %*% 
            BB.time), delta.taylor)
        b.vec.Lambda.lin.star <- sapply(index.vec, survBayes.b.fctn.Lambda.lin, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int.star$b, 
            taylor.time = taylor.coef.time.star$b)
        b.vec.Lambda.quad.a.star <- diag(sapply(index.vec, survBayes.b.fctn.Lambda.quad.a, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int.star$c, 
            taylor.time = taylor.coef.time.star$c))
        b.vec.Lambda.quad.b.star <- cbind(rep(0, KK + 1), diag(sapply(1:KK, 
            survBayes.b.fctn.Lambda.quad.b, time = time, time.int = time.int, 
            time.int.mat = time.int.mat, hr = hr, int = int, 
            int.delta = int.delta, BB.int = BB.int, BB.time = BB.time, 
            KK = KK, taylor.int = taylor.coef.int.star$c, taylor.time = taylor.coef.time.star$c), 
            nrow = KK + 1, ncol = KK))
        b.vec.Lambda.quad.c.star <- cbind(rep(0, KK + 1), rep(0, 
            KK + 1), diag(sapply(1:(KK - 1), survBayes.b.fctn.Lambda.quad.c, 
            time = time, time.int = time.int, time.int.mat = time.int.mat, 
            hr = hr, int = int, int.delta = int.delta, BB.int = BB.int, 
            BB.time = BB.time, KK = KK, taylor.int = taylor.coef.int.star$c, 
            taylor.time = taylor.coef.time.star$c), nrow = KK + 
            1, ncol = KK - 1))
        b.vec.Lambda.quad.d.star <- cbind(rep(0, KK + 1), rep(0, 
            KK + 1), rep(0, KK + 1), diag(sapply(1:(KK - 2), 
            survBayes.b.fctn.Lambda.quad.d, time = time, time.int = time.int, 
            time.int.mat = time.int.mat, hr = hr, int = int, 
            int.delta = int.delta, BB.int = BB.int, BB.time = BB.time, 
            KK = KK, taylor.int = taylor.coef.int.star$c, taylor.time = taylor.coef.time.star$c), 
            nrow = KK + 1, ncol = KK - 2))
        b.vec.Lambda.quad.star <- b.vec.Lambda.quad.a.star + 
            b.vec.Lambda.quad.b.star + t(b.vec.Lambda.quad.b.star) + 
            b.vec.Lambda.quad.c.star + t(b.vec.Lambda.quad.c.star) + 
            b.vec.Lambda.quad.d.star + t(b.vec.Lambda.quad.d.star)
        b.vec.star <- b.vec.cens - b.vec.Lambda.lin.star
        QQ.star <- QQ.hlp + b.vec.Lambda.quad.star
        if (any(eigen(QQ.star, symmetric = TRUE, only.values = TRUE)$values < 
            0)) 
            RRR <- 0
        else {
            LL.star <- t(chol(QQ.star))
#
#               This calculates the R in formula (11) of Rue (2001) for lbh.coef.star
#
            mu.lbh.star <- solve(t(LL.star), solve(LL.star, b.vec.star))
            lbh.int.star <- as.vector(lbh.coef.star %*% BB.int)
            Lambda0.int.star <- int.delta/2 * (exp(lbh.int.star[-(KK - 
                1)]) + exp(lbh.int.star[-1]))
            Lambda0.star <- c(0, cumsum(Lambda0.int.star))
            Lambda0.int.last.star <- (time - int[time.int])/2 * 
                (exp(lbh.int.star[time.int]) + exp(as.vector(lbh.coef.star %*% 
                  BB.time)))
            Lambda0.star <- Lambda0.star[time.int] + Lambda0.int.last.star
            log.pi.comp <- sum(b.vec.cens * (lbh.coef.star - 
                lbh.coef))
            log.pi.comp <- log.pi.comp - sum(hr * (Lambda0.star - 
                Lambda0))
            log.pi.comp <- log.pi.comp - 0.5 * lbh.coef.star %*% 
                QQ.hlp %*% lbh.coef.star
            log.pi.comp <- log.pi.comp + 0.5 * lbh.coef %*% QQ.hlp %*% 
                lbh.coef
            xxx <- (lbh.coef.star - mu.lbh)
            xxx.star <- (lbh.coef - mu.lbh.star)
            log.qq.comp <- 0.5 * (xxx %*% QQ %*% xxx - xxx.star %*% 
                QQ.star %*% xxx.star)
            RRR <- exp(log.pi.comp + log.qq.comp) * prod(diag(LL.star)/diag(LL))
        }
#
#               Now the Metropolis Hastings step will be performed
#
        prop <- ifelse(is.na(RRR), 0, RRR)
        alpha <- min(1, prop)
        u <- runif(1)
        if (u <= alpha) {
            lbh.coef <- lbh.coef.star
            m.h.lbh <- 1
        }
        if (iter > burn.in) 
            lbh.coef.out <- cbind(lbh.coef.out, lbh.coef)
#
#               Perform the sigma.lbh update
#               The update is done for the precisions tau.0=1/sigma.lbh.0
#               and tau.1=1/sigma.lbh.1
#
        lin.form.0 <- 0.5 * (lbh.coef %*% QQ.0 %*% lbh.coef)
        shape.0 <- shape.sigma.lbh.0 + 1/2
        rate.0 <- rate.sigma.lbh.0 + lin.form.0
        sigma.lbh.0 <- 1/rgamma(1, shape = shape.0, rate = 1/rate.0)
        lin.form.1 <- 0.5 * (lbh.coef %*% QQ.ar1 %*% lbh.coef)
        shape.1 <- shape.sigma.lbh.1 + (KK)/2
        rate.1 <- rate.sigma.lbh.1 + lin.form.1
        sigma.lbh.1 <- 1/rgamma(1, shape = shape.1, rate = 1/rate.1)
        if (iter > burn.in) 
            sigma.lbh.out <- cbind(sigma.lbh.out, c(sigma.lbh.0, 
                sigma.lbh.1))
#
#               The following takes care on the random effects
#
        if (!is.null(frailty.values)) {
            lbh.int <- as.vector(lbh.coef %*% BB.int)
            Lambda0.int <- int.delta/2 * (exp(lbh.int[-(KK - 
                1)]) + exp(lbh.int[-1]))
            Lambda0 <- c(0, cumsum(Lambda0.int))
            Lambda0.int.last <- (time - int[time.int])/2 * (exp(lbh.int[time.int]) + 
                exp(as.vector(lbh.coef %*% BB.time)))
            Lambda0 <- Lambda0[time.int] + Lambda0.int.last
            Lambda.x <- exp(X.design %*% beta) * Lambda0
#
#               for lognormal frailty
#
            if (frailty.dist == "gauss") {
                taylor.coef <- survBayes.taylor(alpha.cluster, 
                  delta.taylor)
                b.vec.RE.cens <- as.vector(tapply(cens, frailty.values, 
                  sum))
                b.vec.RE.Lambda <- as.vector(tapply(Lambda.x, 
                  frailty.values, sum))
                b.vec.RE <- b.vec.RE.cens - taylor.coef$b * b.vec.RE.Lambda
                QQ.RE <- rep(1/sigma.RE, n.cl) + taylor.coef$c * 
                  b.vec.RE.Lambda
                mu <- b.vec.RE/QQ.RE
                alpha.cluster.star <- mu + sqrt(1/QQ.RE) * rnorm(n.cl, 
                  0, 1)
#
#               This calculates the R in formula (11) of Rue (2001) for alpha.cluster.star
#
                taylor.coef.star <- survBayes.taylor(alpha.cluster.star, 
                  delta.taylor)
                b.vec.RE.star <- b.vec.RE.cens - taylor.coef.star$b * 
                  b.vec.RE.Lambda
                QQ.RE.star <- rep(1/sigma.RE, n.cl) + taylor.coef.star$c * 
                  b.vec.RE.Lambda
                mu.star <- b.vec.RE.star/QQ.RE.star
                log.pi.comp <- sum(b.vec.RE.cens * (alpha.cluster.star - 
                  alpha.cluster))
                log.pi.comp <- log.pi.comp - sum(b.vec.RE.Lambda * 
                  (exp(alpha.cluster.star) - exp(alpha.cluster)))
                log.pi.comp <- log.pi.comp - 0.5 * sum(alpha.cluster.star^2)/sigma.RE
                log.pi.comp <- log.pi.comp + 0.5 * sum(alpha.cluster^2)/sigma.RE
                xxx <- alpha.cluster.star - mu
                xxx.star <- alpha.cluster - mu.star
                log.qq.comp <- 0.5 * sum(xxx * xxx * QQ.RE) - 
                  0.5 * sum(xxx.star * xxx.star * QQ.RE.star)
                RRR <- exp(log.pi.comp + log.qq.comp) * sqrt(prod(QQ.RE.star/QQ.RE))
#
#               Now the Metropolis Hastings step will be performed
#
                prop <- ifelse(is.na(RRR), 0, RRR)
                alpha <- min(1, prop)
                u <- runif(1)
                if (u <= alpha) {
                  alpha.cluster <- alpha.cluster.star
                  frailty.coef <- alpha.cluster[frailty.values]
                  m.h.alpha <- 1
                }
                if (iter > burn.in) 
                  alpha.cluster.out <- cbind(alpha.cluster.out, 
                    alpha.cluster)
#
#               Perform the sigma.RE update
#
                lin.form.clust <- 0.5 * sum(alpha.cluster * alpha.cluster)
                shape.0 <- shape.sigma.clust + n.cl/2 - 3/2
                rate.0 <- rate.sigma.clust + lin.form.clust
                sigma.RE <- 1/rgamma(1, shape = shape.0, rate = 1/rate.0)
                if (iter > burn.in) 
                  sigma.RE.out <- c(sigma.RE.out, sigma.RE)
            }
#
#               for gamma frailty
#
            else {
                Phi.z <- as.vector(tapply(cens, frailty.values, 
                  sum))
                Psi.z <- as.vector(tapply(Lambda.x, frailty.values, 
                  sum))
                z.cluster <- rgamma(n.cl, shape = (Phi.z + mu.cl), 
                  rate = (Psi.z + mu.cl))
                mm.z <- min(z.cluster[z.cluster > 0])
                z.cluster <- ifelse(z.cluster == 0, mm.z, z.cluster)
                frailty.coef <- log(z.cluster[frailty.values])
                if (iter > burn.in) {
                  z.cluster.out <- cbind(z.cluster.out, z.cluster)
                }
#
#               Perform the mu.cl update
#
                zz <- z.cluster
                AA <- sum(log(zz) - zz)
                rho.star <- prec.tau.cl
#
#               Max of the log-density
#               Using a normal approximation to the cond. posterior at its max as proposal distribution 
#
                newton.raph.max <- 100
                eps.newton.raph <- 1e-05
                tau <- tau.cl
                cnt.newton.raphson <- 0
                while (cnt.newton.raphson < newton.raph.max) {
                  tau.0 <- tau
                  mu.0 <- exp(tau.0)
                  df.0 <- n.cl * mu.0 * (tau.0 + 1 - digamma(mu.0) + 
                    AA/n.cl) - tau.0 * rho.star
                  ddf.0 <- n.cl * mu.0 * (tau.0 + 2 - digamma(mu.0) - 
                    mu.0 * trigamma(mu.0) + AA/n.cl) - rho.star
                  tau <- tau.0 + df.0/(-ddf.0)
                  if (abs(df.0) < eps.newton.raph) 
                    cnt.newton.raphson <- newton.raph.max
                  else cnt.newton.raphson <- cnt.newton.raphson + 
                    1
                }
                mu <- exp(tau)
                ddf <- n.cl * mu * (tau + 2 - digamma(mu) - mu * 
                  trigamma(mu) + AA/n.cl) - rho.star
                sigma.prop <- sqrt((-1)/ddf)
                tau.cl.update <- rnorm(1, tau, sigma.prop)
#
#               Metropolis Hastings Step
#               Make a decision on the proposal
#
                tau.0 <- tau.cl
                tau.1 <- tau.cl.update
                prop.0 <- dnorm(tau.0, tau, sigma.prop)
                prop.1 <- dnorm(tau.1, tau, sigma.prop)
                mu.1 <- exp(tau.1)
                pi.1 <- n.cl * mu.1 * tau.1 - n.cl * log(gamma(mu.1)) + 
                  mu.1 * AA - 0.5 * rho.star * tau.1 * tau.1
                mu.0 <- exp(tau.0)
                pi.0 <- n.cl * mu.0 * tau.0 - n.cl * log(gamma(mu.0)) + 
                  mu.0 * AA - 0.5 * rho.star * tau.0 * tau.0
                RRR <- exp(pi.1 - pi.0) * prop.0/prop.1
                prop <- ifelse(is.na(RRR), 0, RRR)
                alpha <- min(1, prop)
                u <- runif(1)
                if (u <= alpha) {
                  tau.cl <- tau.1
                  mu.cl <- exp(tau.cl)
                  m.h.cl <- 1
                }
                if (iter > burn.in) 
                  mu.cl.out <- c(mu.cl.out, exp(tau.cl))
            }
        }
        if (iter == floor(iter/n.inter) * n.inter) 
            cat(iter, "\n")
        if (is.null(frailty.values)) {
            m.h.performance <- cbind(m.h.performance, c(m.h.beta, 
                m.h.lbh))
        }
        else {
            if (frailty.dist == "gauss") {
                m.h.performance <- cbind(m.h.performance, c(m.h.beta, 
                  m.h.lbh, m.h.alpha))
            }
            else {
                m.h.performance <- cbind(m.h.performance, c(m.h.beta, 
                  m.h.lbh, m.h.cl))
            }
        }
    }
    if (is.null(frailty.values)) {
        res <- list(t.where = int, beta = mcmc(t(beta.out)), 
            lbh.coef = mcmc(t(lbh.coef.out)), sigma.lbh = mcmc(t(sigma.lbh.out)), 
            m.h.performance = apply(m.h.performance, 1, sum))
    }
    else {
        if (frailty.dist == "gauss") {
            res <- list(t.where = int, beta = mcmc(t(beta.out)), 
                lbh.coef = mcmc(t(lbh.coef.out)), sigma.lbh = mcmc(t(sigma.lbh.out)), 
                alpha.cluster = mcmc(t(alpha.cluster.out)), sigma.cluster = mcmc(sigma.RE.out), 
                m.h.performance = apply(m.h.performance, 1, sum))
        }
        else {
            res <- list(t.where = int, beta = mcmc(t(beta.out)), 
                lbh.coef = mcmc(t(lbh.coef.out)), sigma.lbh = mcmc(t(sigma.lbh.out)), 
                z.cluster = mcmc(t(z.cluster.out)), mu.cluster = mcmc(mu.cl.out), 
                m.h.performance = apply(m.h.performance, 1, sum))
        }
    }
    return(res)
}

