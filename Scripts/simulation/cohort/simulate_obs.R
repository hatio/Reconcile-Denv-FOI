#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Simulate age, timings of observations.")

parser$add_argument('-obs', default = '01-processedData/serology/gmt', help = 'Directory of observed data to extract observation years from.')
parser$add_argument('-nogap', default = FALSE, action = "store_true", help = 'Whether to fill in gap years between cohorts.')
parser$add_argument('-ageMin', type = "double", default = 5)
parser$add_argument('-ageMax', type = "double", default = 15)
parser$add_argument('-ageDelta', type = "double", default = 1, help = 'Age difference between each bleed')
parser$add_argument('-num_bleed', type = "integer", default = NULL)
parser$add_argument('-num_bleed_per_indiv', type = "integer", default = 2)


parser$add_argument('-outroot', default = '01-processedData/simcohort/obs', help = 'Root directory to store simulated data.')
parser$add_argument('-dataver', default = 'uniform_5_15', help = 'Name of version of simulated data.')
parser$add_argument('-iFold', type = 'integer', default = 1, help = 'i-th simulation.')
parser$add_argument('-ncpu', type = "integer", default = 1)



# parse arguments
inputArg = parser$parse_args()


set.seed(7350)
runif(inputArg$iFold, 0, 10^6) %>% tail(1) %>% set.seed

outdir = with(inputArg, file.path(outroot, dataver, iFold))
dir.create(outdir, recursive = T)

    # Import observed datasets to get the numbers of individuals under study
    # and timings of their observations

Obs = list.files(inputArg$obs, full.names = T) %>%
    lapply(read_csv, col_types = 'Dcdcdc') %>%
    do.call(what = rbind)

numBleeds = ifelse(
    is.null(inputArg$num_bleed)
    , nrow(Obs)
    , inputArg$num_bleed
)
    

Years = Obs$dateCollect %>%
    format("%Y") %>%
    unique %>%
    as.integer %>%
    sort
if(inputArg$nogap){
    Years = min(Years):max(Years)
}

dat = tibble(year = runif(numBleeds/inputArg$num_bleed_per_indiv, min(Years), max(Years)-1/365 - inputArg$ageDelta)) %>%
    mutate(
        id = seq_along(year)
        , age = runif(n(), inputArg$ageMin, inputArg$ageMax - inputArg$ageDelta)
    ) %>%
    group_by(id) %>%
    nest() %>%
    mutate(data = lapply(data, function(x){
        tibble(
            dateCollect = x$year + seq(0, length.out = inputArg$num_bleed_per_indiv, by = inputArg$ageDelta)
            , age = x$age + seq(0, length.out = inputArg$num_bleed_per_indiv, by = inputArg$ageDelta)
        )
    })) %>%
    unnest(cols = data) %>%
    mutate(
        dateCollect = as.Date(paste0(floor(dateCollect),'-01-01')) + (dateCollect-floor(dateCollect)) * 365
    ) %>%
    select(dateCollect, id, age) %>%
    mutate(
        location = NA
        , titer = NA
        , study = NA
    )

g = ggarrange(
    dat %>%
        mutate(
            year = format(dateCollect, '%Y') %>% as.integer
        ) %>%
        group_by(year) %>%
        summarize(count = n()) %>%
        ggplot(aes(x = year, y = count))+
        geom_col(fill = 'black')
    , dat %>%
        mutate(
            year = format(dateCollect, '%Y') %>% as.integer
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


dat %>%
    write_csv(file.path(outdir, "obs.csv"))
