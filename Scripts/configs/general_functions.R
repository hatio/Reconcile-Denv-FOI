# function to get Wald CI
# : taken from https://towardsdatascience.com/five-confidence-intervals-for-proportions-that-you-should-know-about-7ff5484c024f
waldInterval = function(x, n, conf.level = 0.95){
    p = x/n
    sd = sqrt(p*((1-p)/n))
    z = qnorm(c(
        (1-conf.level) / 2
        , 1 - (1-conf.level)/2
    ))
    ci = p + z*sd
    return(ci)
}

# function to format values to defined number of decimal points
formatValue = function(x, digits = 3){
    x %>%
    round(digits) %>%
    formatC(digits = digits, format = "f", flag = "0")    
}