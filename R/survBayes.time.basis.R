`survBayes.time.basis` <-
function (t, int, KK) 
{
#
#       This function calculates the basis spline functions at the 
#       time points
#
    int.number <- max(which(int < t))
    tm1 <- int[int.number]
    tm2 <- int[int.number - 1]
    tm3 <- int[int.number - 2]
    tp1 <- int[int.number + 1]
    tp2 <- int[int.number + 2]
    tp3 <- int[int.number + 3]
    b4 <- (t - tm1)/(tp3 - tm1) * (t - tm1)/(tp2 - tm1) * (t - 
        tm1)/(tp1 - tm1)
    b3 <- (t - tm2)/(tp2 - tm2) * ((t - tm2)/(tp1 - tm2) * (tp1 - 
        t)/(tp1 - tm1) + (tp2 - t)/(tp2 - tm1) * (t - tm1)/(tp1 - 
        tm1)) + (tp3 - t)/(tp3 - tm1) * (t - tm1)/(tp2 - tm1) * 
        (t - tm1)/(tp1 - tm1)
    b2 <- (t - tm3)/(tp1 - tm3) * (tp1 - t)/(tp1 - tm2) * (tp1 - 
        t)/(tp1 - tm1) + (tp2 - t)/(tp2 - tm2) * ((t - tm2)/(tp1 - 
        tm2) * (tp1 - t)/(tp1 - tm1) + (tp2 - t)/(tp2 - tm1) * 
        (t - tm1)/(tp1 - tm1))
    b1 <- (tp1 - t)/(tp1 - tm3) * (tp1 - t)/(tp1 - tm2) * (tp1 - 
        t)/(tp1 - tm1)
    return(c(rep(0, int.number - 3), b1, b2, b3, b4, rep(0, KK))[1:(KK + 
        1)])
}

