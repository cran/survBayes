"survBayes.chol" <-
function (mat) 
{
#   Choleski-Zerlegung einer Bandmatrix mit Bandweite 1
    DD <- dim(mat)[1]
    OO <- matrix(0, nrow = DD, ncol = DD)
    UU <- matrix(0, nrow = DD, ncol = DD)
    LL <- matrix(0, nrow = DD, ncol = DD)
    UU <- UU + diag(rep(1, DD))
    OO[1, 1] <- mat[1, 1]
    LL[1, 1] <- sqrt(OO[1, 1])
    for (i in 2:DD) {
        UU[i, i - 1] <- mat[i, i - 1]/OO[i - 1, i - 1]
        LL[i, i - 1] <- UU[i, i - 1] * LL[i - 1, i - 1]
        OO[i, i] <- mat[i, i] - mat[i, i - 1] * UU[i, i - 1]
        LL[i, i] <- sqrt(OO[i, i])
        OO[i - 1, i] <- mat[i - 1, i]
    }
    return(LL)
}
