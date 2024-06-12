#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Simulate infection timings of cohort subjects from underlying infection parameters.")

parser$add_argument('-obs', default = '01-processedData/serology/gmt', help = 'Directory of observed data to extract observation timings from.')
parser$add_argument('-params', default = '01-processedData/simcohort/m3/infection_param/mock', help = "File storing the infection parameters (Tau(t), Kappa(a)).")

parser$add_argument('-outroot', default = '01-processedData/simcohort/m3', help = 'Root directory to store simulated data.')
parser$add_argument('-iFold', type = 'integer', default = 1, help = 'i-th simulation.')
parser$add_argument('-ncpu', type = "integer", default = 1)



# parse arguments
inputArg = parser$parse_args()

set.seed(19832)
runif(inputArg$iFold, 0, 10^6) %>% tail(1) %>% set.seed

outdir = with(inputArg, file.path(outroot, iFold))
infectiondir = with(inputArg, file.path(outdir, 'infection'))
dir.create(infectiondir, recursive = T)

    # Import observed datasets to get the numbers of individuals under study
    # and timings of their observations

Obs = list.files(inputArg$obs, pattern = "\\.csv$", full.names = T) %>%
    lapply(read_csv, col_types = 'Dcdcdc')


    # Import KPPH estimates to use as underlying truth

Tau = file.path(inputArg$params, 'Tau.csv') %>% read_csv(col_types = 'dd')
Kappa = file.path(inputArg$params, 'Kappa.csv') %>% read_csv(col_types = 'dd')


    # Simulate to have the same number of individuals as our cohort studies

# functions to draw from truncated distributions
# taken from Nadarajah, 2006
# (https://doi.org/10.18637/jss.v016.c02)
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}




# For each individual in each cohort study
Obs %>%
lapply(function(obs){
    # write headers of infection file for this study
    tibble(
        study = character(0)
        , id = character(0)
        , location = character(0)
        , dateBirth = character(0)
        , infectionAge = character(0)
    ) %>%
    write_csv(file.path(infectiondir, paste0(obs$study[1],'.csv')))
    
    split(obs, obs$id) %>%
        lapply(function(x){
            # birth date
            A.birth = max(x$dateCollect) - max(x$age) * 365.25
            # time lag between start of birth year and birth date
            A.delta = as.numeric(A.birth - as.Date(paste0(format(A.birth, '%Y'), '-01-01')), unit = 'days') / 365.25

            # Risk table
            Risk = Tau %>%
                filter(between(Year
                    , A.birth %>% format('%Y') %>% as.integer
                    , (A.birth + max(Kappa$Age)*365.25) %>% format('%Y') %>% as.integer
                )) %>%
                mutate(Tau = lapply(Tau, function(x){
                    data.frame(
                        Tau = rep(x, each = 2)
                        , weight = c(A.delta, 1 - A.delta)
                    )
                })) %>%
                unnest(cols = Tau) %>%
                tail(-1) %>%
                mutate(Kappa = rep(Kappa$Kappa, each = 2) %>% head(length(Tau))) %>%
                # age at start of age segment
                mutate(ageStart = cumsum(weight) - weight)

            # get age of each infection
            iInfectMax = 4      # max number of infections in a lifetime
            A.infect = sapply(1:iInfectMax, function(i){
                A.infect = rexp(nrow(Risk), with(Risk, Tau * Kappa))
                if(!any(A.infect <= Risk$weight)){ return(NA) }
                # index of age segment which infection occurred
                iInfect = which(A.infect <= Risk$weight)[1]
                # exact age at infection
                Risk$ageStart[iInfect] + A.infect[iInfect]
            }) %>%
            sort

            # write/append infection ages of this individual to file
            x %>%
                head(1) %>%
                select(study, id, location) %>%
                mutate(dateBirth = A.birth) %>%
                mutate(infectAge = paste(A.infect, collapse = ';')) %>%
                write_csv(file.path(infectiondir, paste0(x$study[1],'.csv')), append = T)

        })
    
}) %>%
invisible
