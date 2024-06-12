#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)


source("Scripts/configs/titer_functions.R")


parser = ArgumentParser(description = "Simulate observed titers from infection timings and Ab kinetics parameters.")

parser$add_argument('-obs', default = '01-processedData/serology/gmt', help = 'Directory of observed data to extract observation timings from.')

parser$add_argument('-outroot', default = '01-processedData/simcohort/m3/1', help = 'Root directory to store simulated data.')
parser$add_argument('-ncpu', type = "integer", default = 1)

# Ab kinetics parameters
# : defaults to estimates from Salje, 2018
parser$add_argument('-scenario', default = 'wane', help = 'Name of Ab kinetics scenario.')
parser$add_argument('-Omega', type = "double", default = Omega, nargs = "*", help = "permanent titer")
parser$add_argument('-Gamma', type = "double", default = Gamma, nargs = "*", help = "temporary titer")
parser$add_argument('-Delta', type = "double", default = Delta, nargs = "*", help = "decay rate of temporary titer")
parser$add_argument('-Sigma', type = "double", default = Sigma, help = "assay noise SD in exposed individuals")
parser$add_argument('-Omega0_rel', type = "double", default = 0.2, help = "CXR titer relative to permanent titer")


# parse arguments
inputArg = parser$parse_args()


set.seed(19832)
# runif(inputArg$iFold, 0, 10^6) %>% tail(1) %>% set.seed

outdir = with(inputArg, file.path(outroot, scenario))
dir.create(outdir, recursive = T)

infectiondir = with(inputArg, file.path(outroot, 'infection'))
#dir.create(infectiondir, recursive = T)

    # Import observed datasets to get the numbers of individuals under study
    # and timings of their observations

Obs = list.files(inputArg$obs, pattern = "\\.csv$", full.names = T) %>%
    lapply(read_csv, col_types = 'Dcdcdc')

    # Import infection timings


# Functions to work with these parameters
adjustTiter = function(x){
    log(x/10, 2) + 1
}
adjustTiter.inv = function(x){
    2^x * 10 / 2
}
decayedTiter = function(T.delta, gamma = inputArg$Gamma, omega = inputArg$Omega, delta = inputArg$Delta){
    if(length(omega) > 1){ omega = omega[seq_along(T.delta)] }
    if(length(gamma) > 1){ gamma = gamma[seq_along(T.delta)] }
    if(length(delta) > 1){ delta = delta[seq_along(T.delta)] }
    omega + gamma * exp(-delta * T.delta)
}

    # Simulate observed time points mimicking ones in the cohort studies

Obs %>%
# head(1) %>%
lapply(function(obs){
    
    infections = read_csv(file.path(infectiondir, paste0(obs$study[1],'.csv')), col_types = "cccDc")
    infections = split(infections, infections$id)
    
    split(obs, obs$id) %>%
    # head(3) %>%
    lapply(function(x){
            T.infect = infections[[x$id[1]]] %>%
                with(strsplit(infectionAge, ";")[[1]]) %>%
                as.numeric
            # Based on the timings, generate the true titers using estimates from Salje, 2018
            # true titers at observation points
            Titers = sapply(x$age, function(Time){
                x = decayedTiter(Time - T.infect)
                x[Time < T.infect] = 0
                x = sum(x, na.rm = T)
                ifelse(x > 0, x, with(inputArg, Omega * Omega0_rel))
            })
            
            # generate the observed titers using assay SD estimates from Salje, 2018
            Titers.obs = rnorm(length(Titers), mean = Titers, sd = inputArg$Sigma)

            x %>%
            select(-location) %>%
            mutate(
                Titer.true = Titers %>% adjustTiter.inv
                , titer = Titers.obs %>% adjustTiter.inv
            )
    }) %>%
    do.call(what = rbind) %>%
    write_csv(file.path(outdir, paste0(obs$study[1],'.csv')))

}) %>%
invisible
