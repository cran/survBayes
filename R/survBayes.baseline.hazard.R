"survBayes.baseline.hazard" <-
function (surv.res, type = "log", n.inter = 3, start = NULL, 
    end = NULL, thin = 1) 
{
    int <- surv.res$t.where
    if (!is.null(start)) {
        if (!is.null(end)) 
            lbh.coef <- apply(window(surv.res$lbh.coef, start = start, 
                end = end, thin = thin), 2, mean)
        else lbh.coef <- apply(window(surv.res$lbh.coef, start = start, 
            thin = thin), 2, mean)
    }
    else {
        if (!is.null(end)) 
            lbh.coef <- apply(window(surv.res$lbh.coef, end = end, 
                thin = thin), 2, mean)
        else lbh.coef <- apply(window(surv.res$lbh.coef, thin = thin), 
            2, mean)
    }
    time.max <- max(int)
    KK <- length(int) + 1
    int.left <- int[-(KK - 1)]
    int.right <- int[-1]
    int.delta <- diff(int)
    int.help <- c(int[1] - 2 * int.delta[1], int[1] - int.delta[1], 
        int, time.max + 1, time.max + 2)
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
    if (type != "cum") {
        time <- c(int, time)
        BB <- cbind(BB.int, BB.time)[, order(time)]
        time <- time[order(time)]
        lbh <- as.vector(lbh.coef %*% BB)
        if (type == "log") 
            baseline.hazard <- data.frame(time = time, log.base.haz = lbh)
        if (type == "plain") 
            baseline.hazard <- data.frame(time = time, base.haz = exp(lbh))
    }
    else {
        lbh.int <- as.vector(lbh.coef %*% BB.int)
        Lambda0.int <- int.delta/2 * (exp(lbh.int[-(KK - 1)]) + 
            exp(lbh.int[-1]))
        Lambda0.int <- c(0, cumsum(Lambda0.int))
        Lambda0.time <- NULL
        if (n.inter >= 1) {
            time.int <- sapply(time, function(t, int.left, int.right) {
                which(int <= t & c(int.right, time.max + 1) > 
                  t)
            }, int.left, int.right)
            time.int.mat <- t(sapply(time, function(t, int.right) {
                as.numeric(int.right <= t)
            }, int.right))
            Lambda0.int.last <- (time - int[time.int])/2 * (exp(lbh.int[time.int]) + 
                exp(as.vector(lbh.coef %*% BB.time)))
            Lambda0.time <- Lambda0.int[time.int] + Lambda0.int.last
        }
        time <- c(int, time)
        Lambda0 <- c(Lambda0.int, Lambda0.time)[order(time)]
        time <- time[order(time)]
        baseline.hazard <- data.frame(time = time, cum.base.haz = Lambda0)
    }
    return(baseline.hazard)
}
