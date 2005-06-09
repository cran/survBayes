"survBayes" <-
function (formula = formula(data), data = parent.frame(), burn.in = 1000, 
    number.sample = 1000, max.grid.size = 50, control, control.frailty, 
    seed.set = 100, ...) 
{
    require(survival)
    require(coda)
    call <- match.call()
    m <- match.call(expand = FALSE)
    temp <- c("", "formula", "data")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) 
        terms(formula)
    else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if (NROW(m) == 0) 
        stop("No (non-missing) observations")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, nrow(Y))
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    attr(Terms, "intercept") <- 1
    dropx <- NULL
    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 
        1)
    X <- X[, -1, drop = FALSE]
    type <- attr(Y, "type")
    if (type != "right" && type != "interval") 
        stop("Invalid survival type")
    else {
        if (type == "right") {
            time <- Y[, 1]
            cens <- Y[, 2]
            Z <- cbind(time, cens)
            type.ind <- 1
        }
        if (type == "interval") {
            if (any(Y[, 3] != 3 && Y[, 3] != 0)) 
                stop("Invalid cens status")
            else {
                if (any(Y[, 3] == 3)) 
                  Z <- cbind(Y[, 1:2], Y[, 3])
                else Z <- cbind(Y[, 1], Y[, 3])
            }
            Z[, 2] <- ifelse(Z[, 2] != 1, Z[, 2], NA)
            type.ind <- 2
        }
    }
    frailty.values <- NULL
    frailty.dist <- NULL
    pterms <- sapply(m, inherits, "coxph.penalty")
    if (any(pterms)) {
        pattr <- lapply(m[pterms], attributes)
        temp <- c(attr(Terms, "response"), attr(Terms, "offset"))
        if (length(dropx)) 
            temp <- c(temp, dropx + 1)
        pterms <- pterms[-temp]
        temp <- match((names(pterms))[pterms], attr(Terms, "term.labels"))
        ord <- attr(Terms, "order")[temp]
        if (any(ord > 1)) 
            stop("Penalty terms cannot be in an interaction")
        pcols <- assign[pterms]
        frailty.values <- X[, pcols[[length(pcols)]]]
        X <- X[, -(pcols[[length(pcols)]])]
        frailty.dist <- unlist(strsplit(names(pcols), "\""))[2]
    }
    if (missing(control)) 
        control <- survBayes.control(...)
    if (control$n.inter.miss) 
        control$n.inter <- 100
    if (control$delta.taylor.miss) 
        control$delta.taylor <- 0.1
    if (control$haz.global.miss) 
        control$haz.global <- Inf
    if (control$sigma.lbh.0.miss) 
        control$sigma.lbh.0 <- 100
    if (control$sigma.lbh.1.miss) 
        control$sigma.lbh.1 <- 100
    if (control$prec.beta.miss) 
        control$prec.beta <- 1e-04
    if (control$rate.sigma.lbh.0.miss) 
        control$rate.sigma.lbh.0 <- 1e-04
    if (control$rate.sigma.lbh.1.miss) 
        control$rate.sigma.lbh.1 <- 1e-04
    if (control$shape.sigma.lbh.0.miss) 
        control$shape.sigma.lbh.0 <- 1e-04
    if (control$shape.sigma.lbh.1.miss) 
        control$shape.sigma.lbh.1 <- 1e-04
    if (!is.null(frailty.values)) {
        if (frailty.dist == "gauss") {
            if (missing(control.frailty)) 
                control.frailty <- survBayes.control.lognormal.frailty(...)
            if (control.frailty$sigma.RE.miss) 
                control.frailty$sigma.RE <- 100
            if (control.frailty$rate.sigma.clust.miss) 
                control.frailty$rate.sigma.clust <- 1e-04
            if (control.frailty$shape.sigma.clust.miss) 
                control.frailty$shape.sigma.clust <- 1e-04
        }
        else {
            if (frailty.dist == "gamma") {
                if (missing(control.frailty)) 
                  control.frailty <- survBayes.control.gamma.frailty(...)
                if (control.frailty$prec.tau.cl.miss) 
                  control.frailty$prec.tau.cl <- 1e-04
            }
            else {
                stop("Invalid frailty distribution")
            }
        }
    }
    res <- survBayes.base(int.matrix = Z, type.ind, X.design = X, 
        frailty.values = frailty.values, frailty.dist = frailty.dist, 
        formula = formula, burn.in, number.sample, max.grid.size, 
        data, control, control.frailty, seed.set)
    return(res)
}
