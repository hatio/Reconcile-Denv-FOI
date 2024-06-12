
#   Susceptibility states by age, year
#   .....................
ageMin = 0; ageMax = 30

calcSusceptibility = function(p, i, iMax = 4){
    # i = number of infections acquired
    # p = probability of escaping a particular serotype
    choose(iMax,i) * p^(iMax-i) * (1-p)^i
}

est = extract(fit, c("lambda","binned_kappa_rel"))
len_study_times = (ncol(est$lambda) - T.pre - 1)

    # year, age to estimate susceptibility fractions
    S = expand_grid(
        iTime = 1 : len_study_times
        , age = seq(ageMin, ageMax, by = T.delta)
    )


    # Kappa: iPosterior x age
    if(length(dim(est$binned_kappa_rel)) == 3){
        est$binned_kappa_rel = est$binned_kappa_rel[ , 1, ]
    } else if(is.null(est$binned_kappa_rel)){
        est$binned_kappa_rel = matrix(1, nrow(est$lambda), length(ageBreaks)-1-1)
    }
    est$binned_kappa = cbind(1, est$binned_kappa_rel)
    Kappa = est$binned_kappa[
        , cut(0:ageMax
            , breaks = ageBreaks
            , right = FALSE
            , include.lowest = TRUE
        ) %>% as.integer
    ]
    
    # probability of escaping a particular serotype: (year,age) x iPosterior
    p = 
        # for each individual
        mapply(function(timeBirth, timeObs){
            # for each (year, age) compute presence at risk * kappa
            
            ((Kappa %*%
            # calculate average Kappa of each time interval
            (mapply(function(ageStart, ageEnd){
                    out = encodeTimeInterval(
                        timeBirth + ageStart
                        , min(timeBirth + ageEnd, timeObs)
                        , c(-999, -T.pre:-1, 0:len_study_times) / T.delta
                    )
                    out[out < 0] = 0
                    out
                }
                , ageStart = 0:ageMax
                , ageEnd = (0:ageMax) + 1
                , SIMPLIFY = FALSE
            ) %>%
            do.call(what = rbind))) *
            est$lambda) %>%
            rowSums
        }
        , timeBirth = with(S, iTime - age)
        , timeObs = with(S, iTime)
        , SIMPLIFY = FALSE
    ) %>%
    do.call(what = rbind)
    p = exp(-p)

    # compute probability of being in each susceptibility state
    S = lapply(0:4, function(i){
        Si = calcSusceptibility(p, i) %>%
            apply(1, quantile, c(0.025, 0.5, 0.975), simplify = F) %>%
            do.call(what = rbind) %>%
            as_tibble
        colnames(Si) = c('lowerVal', 'Val', 'upperVal')
        S %>%
            mutate(S = i) %>%
            cbind(Si)
    }) %>%
    do.call(what = rbind)

    # write for later use
    S %>%
        # mutate(year = iTime + 1998 - 1, .before = 'age') %>%
        mutate(year = tail(Times, -T.pre-1)[iTime], .before = 'age') %>%
        write_csv(file.path(outdir, 'S.csv'))
