
#   Susceptibility states by age, year
#   .....................
# ageMin = 3; ageMax = 30

calcSusceptibility = function(p, i, iMax = 4){
    # i = number of infections acquired
    # p = probability of escaping a particular serotype
    choose(iMax,i) * p^(iMax-i) * (1-p)^i
}

est = extract(fit, "lambda")
len_study_times = (ncol(est$lambda) - T.pre - 1)

    # year, age to estimate susceptibility fractions
    S = expand_grid(
        iTime = 1 : len_study_times
        , age = seq(ageMin, ageMax, by = T.delta)
    )

    # for each (year, age) compute presence at risk * kappa
    Kappa = 
        KappaInput %>%
        filter(Age < ageMax + 1) %>%
        with(
            # for each individual
            mapply(function(timeBirth, timeObs){
                Kappa %*%
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
                    , ageStart = Age
                    , ageEnd = Age + 1
                    , SIMPLIFY = FALSE
                ) %>%
                do.call(what = rbind))
            }
            , timeBirth = with(S, iTime - age)
            , timeObs = with(S, iTime)
            , SIMPLIFY = FALSE
        )) %>%
        do.call(what = rbind)
    
        
    # probability of escaping a particular serotype: (year,age) x iPosterior
    p = exp(-Kappa %*% t(est$lambda))

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

