encodeTimeInterval = function(dmin,dmax, Times){
    ((tail(Times, -1) >= dmin) == (head(Times, -1) <= dmax)) * (
        (tail(Times, -1) - dmin) +
        (head(Times, -1) > dmin) * (dmin - head(Times, -1)) +
        (tail(Times, -1) > dmax) * (dmax - tail(Times, -1))
    )
}