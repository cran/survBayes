`survBayes.draw.cond.pw.lin.exp` <-
function (inf.pat, int, Lambda0) 
{
#
#       Patient information inf.pat is a vector of three numbers:
#       t.left,t.right,log.hr
#
    if (is.na(inf.pat[2])) 
        out <- c(time = inf.pat[1], cens = 0)
    else {
        int.left <- int[-length(int)]
        Lambda.x <- exp(inf.pat[3]) * Lambda0
        ll.1 <- approx(int, Lambda.x, inf.pat[1:2])$y
        ll.2 <- 1 - exp(-ll.1)
        ww <- (-1) * log(1 - runif(1, min = ll.2[1], max = ll.2[2]))
        tt <- approx(Lambda.x, int, ww)$y
        out <- c(time = tt, cens = 1)
    }
    return(out)
}

