`survBayes.b.fctn.Lambda.quad.c` <-
function (nth, time, time.int, time.int.mat, hr, int, int.delta, 
    BB.int, BB.time, KK, taylor.int, taylor.time) 
{
#
#       The contribution of relative risk in the nth interval, quadratic part 2. off diagonal
#
    interval.part <- time.int.mat %*% (int.delta/2 * (taylor.int[-(KK - 
        1)] * BB.int[nth, -(KK - 1)] * BB.int[nth + 2, -(KK - 
        1)] + taylor.int[-1] * BB.int[nth, -1] * BB.int[nth + 
        2, -1]))
    rest.part <- (time - int[time.int])/2 * (taylor.int[time.int] * 
        BB.int[nth, time.int] * BB.int[nth + 2, time.int] + taylor.time * 
        BB.time[nth, ] * BB.time[nth + 2, ])
    weighted.Lambda.ind <- (interval.part + rest.part)
    return(sum(hr * weighted.Lambda.ind))
}

