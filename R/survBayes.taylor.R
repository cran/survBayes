"survBayes.taylor" <-
function (u = hh, d = delta) 
{
    e2d <- exp(2 * d)
    edu <- exp(u - d)
    d5 <- d^5
    a <- (e2d - 1) * (d^4 - 15 * u^2 - d^2 * (u * (2 + 5 * u) - 
        5)) + (e2d + 1) * (15 * d * u^2 + d^3 * (2 * u - 5))
    a <- (-1) * a * 3 * edu/(4 * d5)
    b <- 15 * u + d * (d * (1 + d) + 5 * u * (3 + d))
    b <- b + e2d * (d * d * (d - 1) - 5 * u * (3 + d * (d - 3)))
    b <- 3 * edu * b/(2 * d5)
    c <- (3 + d^2) * (e2d - 1) - 3 * d * (1 + e2d)
    c <- c * 15 * edu/(4 * d5)
    return(list(a = a, b = b, c = 2 * c))
}
