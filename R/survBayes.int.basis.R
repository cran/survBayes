`survBayes.int.basis` <-
function (number, int, KK) 
{
#
#       This function calculates the basis spline functions at the 
#       interval points
#
    b1 <- (int[number + 3] - int[number + 2])/(int[number + 3] - 
        int[number]) * (int[number + 3] - int[number + 2])/(int[number + 
        3] - int[number + 1])
    b2 <- (int[number + 2] - int[number])/(int[number + 3] - 
        int[number]) * (int[number + 3] - int[number + 2])/(int[number + 
        3] - int[number + 1]) + (int[number + 4] - int[number + 
        2])/(int[number + 4] - int[number + 1]) * (int[number + 
        2] - int[number + 1])/(int[number + 3] - int[number + 
        1])
    b3 <- (int[number + 2] - int[number + 1])/(int[number + 4] - 
        int[number + 1]) * (int[number + 2] - int[number + 1])/(int[number + 
        3] - int[number + 1])
    return(c(rep(0, number - 1), b1, b2, b3, rep(0, KK - number - 
        1)))
}

