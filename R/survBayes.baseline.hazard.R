`survBayes.baseline.hazard` <-
function (surv.res, type = "log", ci = FALSE, n.inter = 3, start = NULL, 
    end = NULL, thin = 1) 
{
    int <- surv.res$t.where
    time.max <- max(int)
    KK <- length(int) + 1
    int.left <- int[-(KK - 1)]
    int.right <- int[-1]
    int.delta <- diff(int)
    int.help <- c(int[1] - 2 * int.delta[1], int[1] - int.delta[1], 
        int, time.max + int.delta[KK - 2], time.max + int.delta[KK - 
            2])
    BB.int <- sapply(1:(KK - 1), survBayes.int.basis, int.help, 
        KK)
    time <- NULL
    BB.time <- NULL
    if (n.inter >= 1) {
        for (i in 1:n.inter) {
            time <- c(time, int.left + i/(n.inter + 1) * int.delta)
        }
        BB.time <- sapply(time, survBayes.time.basis, int.help, 
            KK)
    }
    if (!is.null(start)) {
        if (!is.null(end)) {
            lbh.coef <- window(surv.res$lbh.coef, start = start, 
                end = end, thin = thin)
        }
        else {
            lbh.coef <- window(surv.res$lbh.coef, start = start, 
                thin = thin)
        }
    }
    else {
        if (!is.null(end)) {
            lbh.coef <- window(surv.res$lbh.coef, end = end, 
                thin = thin)
        }
        else {
            lbh.coef <- window(surv.res$lbh.coef, thin = thin)
        }
    }
    if (type != "cum") {
        time <- c(int, time)
        BB <- cbind(BB.int, BB.time)[, order(time)]
        time <- time[order(time)]
        lbh <- lbh.coef %*% BB
        lbh.mean <- apply(lbh, 2, mean)
        if (ci) 
            lbh.ci <- apply(lbh, 2, quantile, probs = c(0.025, 
                0.975))
        if (type == "log") {
            baseline.hazard <- data.frame(time = time, log.base.haz = lbh.mean)
            if (ci) 
                baseline.hazard <- data.frame(time = time, log.base.haz = lbh.mean, 
                  log.base.haz.lower = lbh.ci[1, ], log.base.haz.upper = lbh.ci[2, 
                    ])
        }
        if (type == "plain") {
            baseline.hazard <- data.frame(time = time, base.haz = exp(lbh.mean))
            if (ci) 
                baseline.hazard <- data.frame(time = time, base.haz = exp(lbh.mean), 
                  base.haz.lower = exp(lbh.ci[1, ]), base.haz.upper = exp(lbh.ci[2, 
                    ]))
        }
    }
    else {
        n.samples <- dim(lbh.coef)[1]
        lbh.int <- lbh.coef %*% BB.int
        Lambda0.int <- (rep(1, n.samples) %*% t(int.delta)) * 
            (exp(lbh.int[, -(KK - 1)]) + exp(lbh.int[, -1]))
        Lambda0.int <- rbind(rep(0, n.samples), apply(Lambda0.int, 
            1, cumsum))
        Lambda0.time <- NULL
        if (n.inter >= 1) {
            warning("n.inter>0 not allowed for type=cum, set to 0")
        }
        time <- int
        Lambda0 <- Lambda0.int
        Lambda0.mean <- apply(Lambda0, 1, mean)
        if (ci) 
            Lambda0.ci <- apply(Lambda0, 1, quantile, probs = c(0.025, 
                0.975))
        baseline.hazard <- data.frame(time = time, cum.base.haz = Lambda0.mean)
        if (ci) 
            baseline.hazard <- data.frame(time = time, cum.base.haz = Lambda0.mean, 
                cum.base.haz.lower = Lambda0.ci[1, ], cum.base.haz.upper = Lambda0.ci[2, 
                  ])
    }
    return(baseline.hazard)
}

