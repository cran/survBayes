"survBayes.base" <-
function (int.matrix, type.ind, X.design, frailty.values, frailty.dist, 
    formula, burn.in, number.sample, max.grid.size, data, control, 
    control.frailty, seed.set) 
{
    set.seed(seed.set)
    n.inter <- control$n.inter                                      # Interval to display number of cycles
    delta.taylor <- control$delta.taylor                            # Width for L2 approximation of exp
#
#
#       Setting up the data
#
#
    if (type.ind == 1) {
        time <- int.matrix[, 1]
        cens <- int.matrix[, 2]
    }
    else {
        time <- ifelse(int.matrix[, 3] == 3, (int.matrix[, 1] + 
            int.matrix[, 2])/2, int.matrix[, 1])
        cens <- int.matrix[, 3]/3
    }
    pp <- ncol(X.design)                                            # number of covariates
    nn <- nrow(X.design)                                            # number of observations
#
#       Setting up the time grid, time intervals
#       The grid is build by an inverse of the KM function of these artificial right censored data. 
#       int.left are the left boundries of the intervals
#       int.right are the right boundries of the intervals
#       int.delta is the length of the single intervals
#       KK is the number of intervals
#
    int <- survBayes.KM.int(time, cens, max.grid.size)
    int.left <- int[-length(int)]
    int.right <- int[-1]
    int.delta <- diff(int)
    KK <- length(int.delta)
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
#
    beta.init <- control$beta.init                                  # Initial values of beta
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
    if (!is.finite(control$haz.global)) {
        haz.global <- sum(cens)/sum(time)
    }
    else {
        haz.global <- control$haz.global
    }
    lbh <- rep(log(haz.global), KK)
    sigma.lbh.0 <- control$sigma.lbh.0
    sigma.lbh.1 <- control$sigma.lbh.1
#
#
#       Parameters for prior distributions needed
#
#
    prec.beta <- control$prec.beta
    Cov.beta.inv <- diag(pp) * prec.beta
    rate.sigma.lbh.0 <- control$rate.sigma.lbh.0
    rate.sigma.lbh.1 <- control$rate.sigma.lbh.1
    shape.sigma.lbh.0 <- control$shape.sigma.lbh.0
    shape.sigma.lbh.1 <- control$shape.sigma.lbh.1
#
#
#       Precision matrix components for random walk
#
#
    QQ.0 <- matrix(0, nrow = KK, ncol = KK)
    QQ.0[1, 1] <- 1
    QQ.1 <- diag(1/int.delta + c(0, 1/int.delta[-KK]))
    for (i in 2:KK) {
        QQ.1[i - 1, i] <- (-1)/int.delta[i - 1]
        QQ.1[i, i - 1] <- (-1)/int.delta[i - 1]
    }
#
#
#       Preparing the sampling
#
#
    iter <- 0
    beta.out <- NULL
    lbh.out <- NULL
    sigma.lbh.out <- NULL
    m.h.performance <- NULL
#
#
#       if frailty term exists, the start must be completed
#
#
    if (is.null(frailty.values)) {
        log.hr <- X.design %*% beta
    }
    else {
        n.cl <- length(unique(frailty.values))                      # number of clusters
#       Inits needed:
#       alpha.cluster contains the cluster specific random effect estimates
#       sigma.RE is the variance of the randomeffect distribution N(0,sigma.re)
        formula.covar <- formula[[3]]
        formula <- as.formula(Surv(time, cens) ~ .)
        formula[[3]] <- formula.covar
        alpha.cluster <- coxph(formula, data)$frail
        sigma.RE <- control.frailty$sigma.RE
#       Parameters for prior distributions needed
        rate.sigma.clust <- control.frailty$rate.sigma.clust
        shape.sigma.clust <- control.frailty$shape.sigma.clust
#       Preparing the sampling
        alpha.cluster.out <- NULL
        sigma.RE.out <- NULL
        log.hr <- alpha.cluster[frailty.values] + X.design %*% 
            beta
    }
#
#
#       Starting the sampling
#
#
    sample.total <- burn.in + number.sample
    while (iter < sample.total) {
        iter <- iter + 1
        m.h.beta <- 0
        m.h.lbh <- 0
        m.h.alpha <- 0
#
#       make the data augmentation
#
        if (type.ind == 2) {
            if (is.null(frailty.values)) {
                log.hr <- X.design %*% beta
            }
            else {
                log.hr <- alpha.cluster[frailty.values] + X.design %*% 
                  beta
            }
            inf.pat.mat <- cbind(int.matrix[, 1:2], log.hr)
            time.aug <- apply(inf.pat.mat, 1, survBayes.draw.cond.pw.lin.exp, 
                lbh = lbh, int = int)
            time <- time.aug[1, ]
            cens <- time.aug[2, ]
        }
#
#
#               Perform the beta update
#
#
        log.Lambda.0 <- log(apply(matrix(time, ncol = 1), 1, 
            survBayes.Lambda0, left = int.left, right = int.right, 
            ln.lambda0 = lbh))
        if (is.null(frailty.values)) {
            coef.update <- survBayes.poisson.update(beta, X.design, 
                cens, log.Lambda.0, Cov.beta.inv)
        }
        else {
            coef.update <- survBayes.poisson.update(beta, X.design, 
                cens, alpha.cluster[frailty.values] + log.Lambda.0, 
                Cov.beta.inv)
        }
        u <- runif(1)
        if (u <= coef.update$alpha) {
            beta <- coef.update$beta.star
            m.h.beta <- 1
        }
        if (iter > burn.in) 
            beta.out <- cbind(beta.out, beta)
#
#
#               Perform the lbh update following section 3.1.2 of Rue (2001)
#
#
        ebx <- exp(X.design %*% beta)
        taylor.coef <- survBayes.taylor(lbh, delta.taylor)
#
#
#               The b-term in fomula (4) of Rue 2001 is calculated
#
#
        b.vec.cens <- sapply(1:KK, survBayes.numb.events.int, 
            time = time, cens = cens, left = int.left, right = int.right)
        b.vec.Lambda <- sapply(1:KK, survBayes.b.fctn.Lambda, 
            time = time, ebx = ebx, left = int.left, right = int.right)
        b.vec <- b.vec.cens - taylor.coef$b * b.vec.Lambda
#
#               The Q-term in fomula (4) of Rue (2001) is calculated
#
        QQ.hlp <- QQ.0/sigma.lbh.0 + QQ.1/sigma.lbh.1
        QQ <- QQ.hlp + diag(taylor.coef$c * b.vec.Lambda)
#
#               Sampling a new candidate for lbh following section 3.1.2 of Rue (2001)
#
        LL <- survBayes.chol(QQ)
        mu <- solve(t(LL), solve(LL, b.vec))
        lbh.star <- mu + solve(t(LL), rnorm(KK))
#
#               Metropolis Hastings step to accept the new candidate
#
        taylor.coef.star <- survBayes.taylor(lbh.star, delta.taylor)
        b.vec.star <- b.vec.cens - taylor.coef.star$b * b.vec.Lambda
        QQ.star <- QQ.hlp + diag(taylor.coef.star$c * b.vec.Lambda)
        LL.star <- survBayes.chol(QQ.star)
        mu.star <- solve(t(LL.star), solve(LL.star, b.vec.star))
#
#
#               This calculates the R in formula (11) of Rue (2001) for lbh.star
#
#
        log.pi.comp <- sum(b.vec.cens * (lbh.star - lbh))
        log.pi.comp <- log.pi.comp - sum(b.vec.Lambda * (exp(lbh.star) - 
            exp(lbh)))
        log.pi.comp <- log.pi.comp - 0.5 * lbh.star %*% QQ.hlp %*% 
            lbh.star
        log.pi.comp <- log.pi.comp + 0.5 * lbh %*% QQ.hlp %*% 
            lbh
        xxx <- (lbh.star - mu)
        xxx.star <- (lbh - mu.star)
        log.qq.comp <- 0.5 * (xxx %*% QQ %*% xxx - xxx.star %*% 
            QQ.star %*% xxx.star)
        RRR <- exp(log.pi.comp + log.qq.comp) * prod(diag(LL.star)/diag(LL))
#
#
#               Now the Metropolis Hastings step will be performed
#
#
        prop <- ifelse(is.na(RRR), 0, RRR)
        alpha <- min(1, prop)
        u <- runif(1)
        if (u <= alpha) {
            lbh <- lbh.star
            m.h.lbh <- 1
        }
        if (iter > burn.in) 
            lbh.out <- cbind(lbh.out, lbh)
#
#
#               Perform the sigma.lbh update
#               The update is done for the precisions tau.0=1/sigma.lbh.0
#               and tau.1=1/sigma.lbh.1
#
#
        lin.form.0 <- 0.5 * (lbh %*% QQ.0 %*% lbh)
        shape.0 <- shape.sigma.lbh.0 + 1/2
        rate.0 <- rate.sigma.lbh.0 + lin.form.0
        sigma.lbh.0 <- 1/rgamma(1, shape = shape.0, rate = rate.0)
        lin.form.1 <- 0.5 * (lbh %*% QQ.1 %*% lbh)
        shape.1 <- shape.sigma.lbh.1 + (KK - 1)/2
        rate.1 <- (rate.sigma.lbh.1 + lin.form.1)
        sigma.lbh.1 <- 1/rgamma(1, shape = shape.1, rate = rate.1)
        if (iter > burn.in) 
            sigma.lbh.out <- cbind(sigma.lbh.out, c(sigma.lbh.0, 
                sigma.lbh.1))
#
#
#               The following takes care on the random effects using a Rue (2001)
#               and a MH update
#
#
#
        if (!is.null(frailty.values)) {
            taylor.coef <- survBayes.taylor(alpha.cluster, delta.taylor)
            Lambda.0 <- apply(matrix(time, ncol = 1), 1, survBayes.Lambda0, 
                left = int.left, right = int.right, ln.lambda0 = lbh)
            Lambda.x <- exp(X.design %*% beta) * Lambda.0
            b.vec.RE.cens <- as.vector(tapply(cens, frailty.values, 
                sum))
            b.vec.RE.Lambda <- as.vector(tapply(Lambda.x, frailty.values, 
                sum))
            b.vec.RE <- b.vec.RE.cens - taylor.coef$b * b.vec.RE.Lambda
            QQ.RE <- rep(1/sigma.RE, n.cl) + taylor.coef$c * 
                b.vec.RE.Lambda
            mu <- b.vec.RE/QQ.RE
            alpha.cluster.star <- mu + sqrt(1/QQ.RE) * rnorm(n.cl, 
                0, 1)
#
#
#               This calculates the R in formula (11) of Rue (2001) for alpha.cluster.star
#
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
            log.qq.comp <- 0.5 * sum(xxx * xxx * QQ.RE) - 0.5 * 
                sum(xxx.star * xxx.star * QQ.RE.star)
            RRR <- exp(log.pi.comp + log.qq.comp) * sqrt(prod(QQ.RE.star/QQ.RE))
#
#
#               Now the Metropolis Hastings step will be performed
#
#
            prop <- ifelse(is.na(RRR), 0, RRR)
            alpha <- min(1, prop)
            u <- runif(1)
            if (u <= alpha) {
                alpha.cluster <- alpha.cluster.star
                m.h.alpha <- 1
            }
            if (iter > burn.in) 
                alpha.cluster.out <- cbind(alpha.cluster.out, 
                  alpha.cluster)
#
#
#               Perform the sigma.RE update
#
#
            lin.form.clust <- 0.5 * sum(alpha.cluster * alpha.cluster)
            shape.0 <- shape.sigma.clust + n.cl/2
            rate.0 <- rate.sigma.clust + lin.form.0
            sigma.RE <- 1/rgamma(1, shape = shape.0, rate = rate.0)
            if (iter > burn.in) 
                sigma.RE.out <- c(sigma.RE.out, sigma.RE)
        }
#
#
#
#
        if (iter == floor(iter/n.inter) * n.inter) 
            cat(iter, "\n")
        if (is.null(frailty.values)) {
            m.h.performance <- cbind(m.h.performance, c(m.h.beta, 
                m.h.lbh))
        }
        else {
            m.h.performance <- cbind(m.h.performance, c(m.h.beta, 
                m.h.lbh, m.h.alpha))
        }
    }
    if (is.null(frailty.values)) {
        res <- list(t.where = int.left, lbh = mcmc(t(lbh.out)), 
            beta = mcmc(t(beta.out)), sigma.lbh = mcmc(t(sigma.lbh.out)), 
            m.h.performance = apply(m.h.performance, 1, sum))
    }
    else {
        res <- list(t.where = int.left, lbh = mcmc(t(lbh.out)), 
            beta = mcmc(t(beta.out)), sigma.lbh = mcmc(t(sigma.lbh.out)), 
            alpha.cluster = mcmc(t(alpha.cluster.out)), sigma.cluster = mcmc(sigma.RE.out), 
            m.h.performance = apply(m.h.performance, 1, sum))
    }
    return(res)
}
