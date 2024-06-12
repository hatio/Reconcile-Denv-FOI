dat = 
    list.files(indir.sero, full.names = T) %>%
    lapply(read_csv) %>%
    do.call(what = rbind) %>%
    filter(
        !is.na(titer)
        , between(age, ageMin.sero, ageMax.sero)
    ) %>%
    mutate(
        # convert dateCollect to fractional years
        year = as.numeric(format(dateCollect,'%Y')) +
            as.numeric(dateCollect - as.Date(paste0((format(dateCollect,'%Y')),'-01-01')))/365.25
    )

# seroincidence data
dat.inc = dat %>%
    select(study, id, age, year, titer) %>%
    group_by(study, id) %>%
    nest() %>%
    mutate(data = lapply(data, function(x){
        x = with(x, tibble(
            age1 = head(age, -1)
            , year1 = head(year, -1)
            , year2 = tail(year, -1)
            , titer1 = head(titer, -1)
            , titer2 = tail(titer, -1)
        )) %>%
        filter(titer1 < max(titer.threshold))
    })) %>%
    group_by(study) %>%
    mutate(idOrder = seq_along(id)) %>%
    unnest(cols = data) %>%
    mutate(year0 = year1 - age1) %>%
    ungroup

# seroprevalence data: include only one bleed per individual
if(!is.na(iFold)){
    # include one random draw per individual
    set.seed(9472224 * 65)
    set.seed(runif(iFold)[iFold] * 10^7)
    dat.prev = dat %>%
        group_by(study, id) %>%
        filter(sample(seq_along(id))==1)

    # sample only one individual per clusterId in KFCS
    if('kfcs' %in% dat.prev$study){
        set.seed(4268877 * ifelse(is.na(iFold), 98, iFold * 43))
        dat.prev = split(dat.prev, dat.prev$study %in% c('kfcs'))
        dat.prev = dat.prev[['TRUE']] %>%
            mutate(clusterId = str_extract(id, "[0-9]{3}KF")) %>%
            group_by(clusterId) %>%
            filter(rank(runif(n())) == 1) %>%
            ungroup %>%
            select(-clusterId) %>%
            rbind(dat.prev[['FALSE']])        
    }
} else {
    dat.prev = dat
}

# determine seropositivity
for (iThresh in 1:length(titer.threshold)){
    # seroprevalence
    dat.prev[[paste0('posAlt.', iThresh)]] = dat.prev$titer >= titer.threshold[iThresh]
    # seroincidence
    dat.inc[[paste0('pos1Alt.', iThresh)]] = dat.inc $titer1>= titer.threshold[iThresh]
    dat.inc[[paste0('pos2Alt.', iThresh)]] = dat.inc $titer2>= titer.threshold[iThresh]
}


# times of interest for joint model
Times = mapply(function(t1,t2){
        seq(floor(t1), ceiling(t2), by = T.delta)
    }
    , t1 = dat.inc$year1 - T.pre*T.delta
    , t2 = dat.inc$year2
    , SIMPLIFY = FALSE
) %>%
Reduce(f = union)
Times = c( with(dat.inc, floor(min(year1 - age1))), Times) %>%
    unique %>%
    sort

# write for later use
as.character(Times) %>%
    writeLines(file.path(outdir, 'Times.txt'))

