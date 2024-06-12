
dat = 
    list.files(indir, full.names = T) %>%
    lapply(read_csv) %>%
    do.call(what = rbind)

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
    select(study, id, age, year, pos)
    
# include only one bleed per individual
if(!is.na(iFold)){
    # include one random draw per individual
    set.seed(9472224 * 65)
    set.seed(runif(iFold)[iFold] * 10^7)
    dat = dat %>%
        group_by(study, id) %>%
        filter(sample(seq_along(id))==1)

    # sample only one individual per clusterId in KFCS
    if('kfcs' %in% dat$study){
        set.seed(4268877 * ifelse(is.na(iFold), 98, iFold * 43))
        dat = split(dat, dat$study %in% c('kfcs'))
        dat = dat[['TRUE']] %>%
            mutate(clusterId = str_extract(id, "[0-9]{3}KF")) %>%
            group_by(clusterId) %>%
            filter(rank(runif(n())) == 1) %>%
            ungroup %>%
            select(-clusterId) %>%
            rbind(dat[['FALSE']])        
    }
}

Studies = c("kps1", "kps2", "kps3", "kfcs")
Studies = unique(c(Studies, dat$study))
print(Studies)
dat =
    dat %>%
    mutate(study = factor(study, levels = Studies)) %>%
    arrange(study)


g = ggarrange(
    dat %>%
        mutate(
            year = floor(year)
        ) %>%
        group_by(year) %>%
        summarize(count = n()) %>%
        ggplot(aes(x = year, y = count))+
        geom_col(fill = 'black')
    , dat %>%
        mutate(
            year = floor(year)
        ) %>%
        group_by(age, year) %>%
        summarize(count = n()) %>%
        ggplot(aes(x = year, y = age, group = year))+
        geom_violin(fill = 'black')
    , ncol = 1
    , nrow = 2
    , heights = c(1,2)
    , align = 'v'
)
ggsave(g
    , filename = file.path(outdir, "data.pdf")
    , width = 6
    , height = 4
)

# times of interest
Times = seq(floor(min(dat$year)) - T.pre*T.delta, ceiling(max(dat$year)), by = T.delta)
Times = c( with(dat, floor(min(year - age))), Times)

# write for later use
as.character(Times) %>%
    writeLines(file.path(outdir, 'Times.txt'))
