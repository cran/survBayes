"survBayes.numb.events.int" <-
function (nth, time, cens, left, right) 
{
#
#       Calculates the number of events in the nth grid interval
#
    left.ind <- time >= left[nth]
    right.ind <- time < right[nth]
    ind <- left.ind & right.ind
    return(sum(cens * ind))
}
