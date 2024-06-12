
if(ppos.ver == "set"){
    # pre-specified P(pos)
    p_pos = cbind(
        c( ppos0, ppos1.long, 1)
        , matrix(c( ppos0, ppos1.short, 1), 3, length(Studies))
    )

# inference on simulated data
} else {

    Sigma = ifelse(grepl('noise', ppos.ver), Sigma, 0)
    Omega0 = 0.2 * ifelse(grepl('noise', ppos.ver) & !grepl('noisepos', ppos.ver), Omega, 0)
    Delta = ifelse(grepl('wane', ppos.ver), Delta, 0)
    Omega_long  = Omega + Gamma * exp(-Delta * 100)

    # long-term behavior (seroprevalence)
    p_pos = matrix(c(
          1 - pnorm(adjustTiter(titer.threshold), Omega0    , Sigma)        # naive
        , 1 - pnorm(adjustTiter(titer.threshold), Omega_long, Sigma)        # monotypic
        , 1                                                                 # multitypic
    ), 3, 1)

    # if fitting to seroincidence data
    if('year2' %in% names(dat)){
        
        deltaT =
            split(dat, dat$study) %>%
            sapply(function(x){
                if(nrow(x) == 0){ return(0)}
                with(x, mean(year2 - year1)/2)
            })

        Omega_short = Omega + Gamma * exp(-Delta * deltaT)

        p_pos = matrix(p_pos, 3, 1 + length(Studies))
        p_pos[2, -1] = 1 - pnorm(adjustTiter(titer.threshold), Omega_short, Sigma)
    }
}
