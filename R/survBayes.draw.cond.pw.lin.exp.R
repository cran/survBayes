"survBayes.draw.cond.pw.lin.exp" <-
function (inf.pat, lbh, int) 
{
#
#
#       Patient information inf.pat is a vector of three numbers:
#       t.left,t.right,log.hr
#
#
    if (is.na(inf.pat[2])) 
        out <- c(time = inf.pat[1], cens = 0)
    else {
        int.left <- int[-length(int)]
        int.right <- int[-1]
        int.delta <- int.right - int.left
#
#       obs.t is in the last interval [int[KK],Inf[
#
        if (inf.pat[2] <= max(int.left)) {
            inf.pat[2] <- min(c(inf.pat[2], max(int.left)))
            Lambda.x <- c(0, cumsum(exp(inf.pat[3] + lbh) * int.delta))
            ll.1 <- approx(int, Lambda.x, inf.pat[1:2])$y
            ll.2 <- 1 - exp(-ll.1)
            ww <- (-1) * log(1 - runif(1, min = ll.2[1], max = ll.2[2]))
            tt <- approx(Lambda.x[Lambda.x < Inf], int.left, 
                ww)$y
            out <- c(time = tt, cens = 1)
        }
        if (inf.pat[1] >= max(int.left)) {
            out <- c(time = 0.5 * (inf.pat[1] + inf.pat[2]), 
                cens = 1)
        }
        if (inf.pat[1] < max(int.left) & inf.pat[2] >= max(int.left)) {
            inf.pat[2] <- min(c(inf.pat[2], max(int.left)))
            Lambda.x <- c(0, cumsum(exp(inf.pat[3] + lbh) * int.delta))
            ll.1 <- approx(int, Lambda.x, inf.pat[1:2])$y
            ll.2 <- 1 - exp(-ll.1)
            ww <- (-1) * log(1 - runif(1, min = ll.2[1], max = ll.2[2]))
            tt <- approx(Lambda.x[Lambda.x < Inf], int.left, 
                ww)$y
            out <- c(time = tt, cens = 1)
        }
    }
    return(out)
}
