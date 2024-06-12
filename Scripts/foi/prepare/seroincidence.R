

dat = 
    list.files(indir, full.names = T) %>%
    lapply(read_csv) %>%
    do.call(what = rbind) %>%
    mutate(study = ifelse(is.na(study), "NA", study)) %>%
    filter(!(study %in% studyExclude))

dat = dat %>%
    filter(
        !is.na(titer)
        , between(age, ageMin, ageMax)
    ) %>%
    mutate(pos = titer >= titer.threshold) %>%
    mutate(
        # convert dateCollect to fractional years
        year = as.numeric(format(dateCollect,'%Y')) +
            as.numeric(dateCollect - as.Date(paste0((format(dateCollect,'%Y')),'-01-01')))/365.25
    ) %>%
    select(study, id, age, year, pos) %>%
    group_by(study, id) %>%
    nest() %>%
    mutate(data = lapply(data, function(x){
        x = with(x, tibble(
            age1 = head(age, -1)
            , year1 = head(year, -1)
            , year2 = tail(year, -1)
            , pos1 = head(pos, -1)
            , pos2 = tail(pos, -1)
        )) %>%
        filter(!pos1)
    })) %>%
    group_by(study) %>%
    mutate(idOrder = seq_along(id)) %>%
    unnest(cols = data)
    
dat = dat %>%
    mutate(year0 = year1 - age1)

Studies = c("kps1", "kps2", "kps3", "kfcs")
Studies = unique(c(Studies, dat$study))
print(Studies)
dat =
    dat %>%
    mutate(study = factor(study, levels = Studies)) %>%
    arrange(study)


# times covered by the data
Times.data =
    with(dat, c(year1,year2)) %>%
    exec(.fn = function(x){
        (
            floor(x) +
            floor((x - floor(x))/T.delta) * T.delta
        ) %>%
        unique %>%
        sort  
    })

# times of interest
Times = mapply(function(t1,t2){
        seq(floor(t1), ceiling(t2), by = T.delta)
    }
    , t1 = dat$year1 - T.pre*T.delta
    , t2 = dat$year2
    , SIMPLIFY = FALSE
) %>%
Reduce(f = union) %>%
sort
Times = c( with(dat, floor(min(year1 - age1))), Times)

# write for later use
as.character(Times.data) %>%
    writeLines(file.path(outdir, 'Times_data.txt'))

as.character(Times) %>%
    writeLines(file.path(outdir, 'Times.txt'))

