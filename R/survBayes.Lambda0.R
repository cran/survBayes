"survBayes.Lambda0" <-
function (time, left, right, ln.lambda0) 
{
    left.min <- pmin(left, time)
    right.min <- pmin(right, time)
    min.diff <- right.min - left.min
    Lambda0.int <- exp(ln.lambda0) * min.diff
    sum(Lambda0.int)
}
