
#   FOI
#   ...

fitsum = summary(fit)$summary
iRow = grepl('^lambda\\[', rownames(fitsum))

est = data.frame(
        year = head(Times,-1)
        , val = fitsum[iRow,'50%']
        , upperVal = fitsum[iRow,'97.5%']
        , lowerVal = fitsum[iRow,'2.5%']
    )
write_csv(est, file.path(outdir, 'foi.csv'))

g = 
    est %>%
    rbind(tail(est,1) %>% mutate(year = max(est$year+1))) %>%
    ggplot(aes(x = year, y = val))+
    geom_step()+
    geom_step(aes(y = upperVal), linetype = 3)+
    geom_step(aes(y = lowerVal), linetype = 3)+
    ggtitle(inputArg$stanfile)+
    ylab("Per-serotype FOI")+
    coord_cartesian(xlim = c(Times[2], max(Times)))
ggsave(g
    , filename = file.path(outdir, 'foi.pdf')
    , width = 6
    , height = 4
)



#   other parameters
#   ................

est = data.frame(
        param = rownames(fitsum)[!iRow]
        , val = fitsum[!iRow,'50%']
        , upperVal = fitsum[!iRow,'97.5%']
        , lowerVal = fitsum[!iRow,'2.5%']
        , n_eff = fitsum[!iRow,'n_eff']
        , Rhat = fitsum[!iRow,'Rhat']
    )
write_csv(est, file.path(outdir, 'est.csv'))
