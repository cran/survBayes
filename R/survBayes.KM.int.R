"survBayes.KM.int" <-
function (time, cens, n.gr) 
{
#
#       This function initializes the time grid over which the hazard function
#       is constant based on the inverse image of an equally spaced grid on [0,1]
#       using the inverse function of the KM estimate (linearly interpolated)
#
    sv <- survfit(Surv(time, cens) ~ 1)
    x.sv <- c(0, sv$time[sv$n.event > 0])
    y.sv <- c(1, sv$surv[sv$n.event > 0])
    y.grid <- seq(min(y.sv), 1, length = n.gr)
    t.grid <- approx(y.sv, x.sv, y.grid)$y
    return(t.grid[order(t.grid)])
}
