#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Compute susceptibility fractions from underlying infection parameters.")

parser$add_argument('params', help = "File storing the infection parameters (Tau(t), Kappa(a)).")

parser$add_argument('-ageMin', type = "double", default = 3)
parser$add_argument('-ageMax', type = "double", default = 30)
parser$add_argument('-T.delta', type = "double", default = 1)
parser$add_argument('-yearStart', type = "double", default = 1998)


# parse arguments
inputArg = parser$parse_args()


calcSusceptibility = function(p, i, iMax = 4){
    # i = number of infections acquired
    # p = probability of escaping a particular serotype
    choose(iMax,i) * p^(iMax-i) * (1-p)^i
}



#   Susceptibility states by age, year
#   .....................

Tau = file.path(inputArg$params, 'Tau.csv') %>% read_csv
Kappa = file.path(inputArg$params, 'Kappa.csv') %>% read_csv


S = expand_grid(
    # year, age to estimate susceptibility fractions
    Year = inputArg$yearStart : max(Tau$Year)
    , Age = with(inputArg, seq(ageMin, ageMax, by = T.delta))
) %>%
mutate(Si = mapply(function(y,a){
        # probability of escaping a particular serotype
        p = Tau %>%
            filter(between(Year, y-a, y)) %>%
            cbind(Kappa %>% filter(Age <= a)) %>%
            with(exp(-sum(Tau * Kappa)))
        # compute probability of being in each susceptibility state
        tibble(S = 0:4) %>%
            mutate(Val = sapply(S, calcSusceptibility, p = p))
    }
    , y = Year
    , a = Age
    , SIMPLIFY = FALSE
)) %>%
unnest(cols = Si)

# write for later use
S %>%
    rename(year = Year, age = Age) %>%
    mutate(lowerVal = Val, upperVal = Val) %>%
    write_csv(file.path(inputArg$params, 'S.csv'))
