"survBayes.b.fctn.Lambda" <-
function (nth, time, ebx, left, right) 
{
#
#       The contribution of relative risk in the nth interval
#
#
    left.min <- pmin(left[nth], time)
    right.min <- pmin(right[nth], time)
    min.diff <- right.min - left.min
    Lambda.int <- ebx * min.diff
    b.Lambda <- sum(Lambda.int)
    return(b.Lambda)
}
