"survBayes.taylor" <-
function (l0, d) 
{
    ep <- exp(l0 + d)
    em <- exp(l0 - d)
    d5 <- d^5
    a <- (5 * l0^2 * (d^2 + 3 * d + 3) + 2 * l0 * d^2 * (d + 
        1) - d^2 * (d^2 + 5 * d + 5)) * em + (-5 * l0^2 * (d^2 - 
        3 * d + 3) + 2 * l0 * d^2 * (d - 1) + d^2 * (d^2 - 5 * 
        d + 5)) * ep
    a <- -a * 3/(4 * d5)
    b <- (5 * l0 * (d^2 + 3 * d + 3) + d^2 * (d + 1)) * em + 
        (-5 * l0 * (d^2 - 3 * d + 3) + d^2 * (d - 1)) * ep
    b <- b * 3/(2 * d5)
    c <- (d^2 + 3 * d + 3) * em - (d^2 - 3 * d + 3) * ep
    c <- -c * 15/(2 * d5)
    return(list(a = a, b = b, c = c))
}
