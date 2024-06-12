#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(rstan)

parser = ArgumentParser(description = "Fit FOI from seroconversion data.")

parser$add_argument('indir')
parser$add_argument('-outbase', required = T)
parser$add_argument('-stanver', default = "time")
parser$add_argument('-studyExclude', default = "")

parser$add_argument('-titer.threshold', type = "double", default = 20)
parser$add_argument('-model.ver', default = "simple")
parser$add_argument('-ppos0', type = "double", default = 0, help = "Test positive probability at S=0")
parser$add_argument('-ppos1.short', type = "double", default = 1, help = "Short-term test positive probability at S=1")
parser$add_argument('-ppos1.long' , type = "double", default = 1, help =  "Long-term test positive probability at S=1")

parser$add_argument('-ageMin', type = "double", default = 3)
parser$add_argument('-ageMax', type = "double", default = 30)
parser$add_argument('-T.delta', type = "double", default = 1)
parser$add_argument('-T.pre', type = "integer", default = 10, help = "Number of time-varying FOI pre-study period before assuming constant FOI.")
parser$add_argument('-kappafile', default = "01-processedData/simcohort/m3/infection_param/mock/Kappa.csv", help = "File containing Kappa(a) values to give to the model.")

parser$add_argument('-ncpu', type = "integer", default = 3)
parser$add_argument('-thin_factor', type = "integer", default = 5)
parser$add_argument('-n_iter', type = "integer", default = 2000)
parser$add_argument('-n_warmup', type = "integer", default = 200)
parser$add_argument('-overwriteRDS', action = "store_true", default = FALSE)


# parse arguments
inputArg = parser$parse_args()
attach(inputArg)


stanfile = file.path("Scripts/foi/seroincidence", stanver,"model.stan")
message(paste("Using stanfile:", stanfile))


outdir = file.path(outbase, stanver, model.ver)
if(grepl('^set', model.ver)){
    outdir = file.path(outdir, paste(c(ppos0, ppos1.long, ppos1.short), collapse = '_'))
}
dir.create(outdir, recursive = T)
stanRDS = file.path(outdir, "fit.RDS")

# load necessary functions
source("Scripts/configs/encoding_functions.R")
source("Scripts/configs/titer_functions.R")



    #   Data
    #   ....

source("Scripts/foi/prepare/seroincidence.R")


    #   Fit model
    #   .........

KappaInput =
    read_csv(kappafile) %>%
    filter(Age < ageMax + 1) %>%
    mutate(Kappa = ifelse( rep(grepl('kappa', model.ver), length(Kappa))
        , Kappa
        , 1
    ))

if(!file.exists(stanRDS) | overwriteRDS){
    
    # presence matrices, for each interval (interval x time)
    # encoding Kappa * time spent in the time interval
    Kappa = 
        KappaInput %>%
        with(
            # for each individual
            mapply(function(timeBirth, timeStart, timeEnd){
                Kappa %*%
                # calculate average Kappa of each time interval
                (mapply(function(ageStart, ageEnd){
                        out = encodeTimeInterval(
                            max(timeBirth + ageStart, timeStart)
                            , min(timeBirth + ageEnd, timeEnd)
                            , Times
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
            , timeBirth = with(dat, year1 - age1)
            , timeStart = with(dat, year1)
            , timeEnd = with(dat, year2)
            , SIMPLIFY = FALSE
        )) %>%
        do.call(what = rbind)

    preKappa = 
        KappaInput %>%
        filter(Age < ageMax + 1) %>%
        with(
            # for each individual
            mapply(function(timeBirth, timeStart, timeEnd){
                Kappa %*%
                # calculate average Kappa of each time interval
                (mapply(function(ageStart, ageEnd){
                        out = encodeTimeInterval(
                            max(timeBirth + ageStart, timeStart)
                            , min(timeBirth + ageEnd, timeEnd)
                            , Times
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
            , timeBirth = with(dat, year0)
            , timeStart = with(dat, year0)
            , timeEnd = with(dat, year1)
            , SIMPLIFY = FALSE
        )) %>%
        do.call(what = rbind)
    
    # Probability of testing positive
    ppos.ver = gsub('_*kappa[A-Za-z]+', '', model.ver)
    source("Scripts/foi/prepare/ppos.R")

    # form input for stan
    stanInput = list(
        num_times = length(Times) - 1
        , num_intervals = nrow(dat)
        , num_study = length(Studies)
        , num_intervals_study = factor(dat$study, levels = Studies) %>% table
        , obs = as.integer(dat$pos2)
        , Presence = Kappa
        , prePresence = preKappa
        , p_pos = p_pos
    )

    set.seed(ncpu*200)
    message(paste("Initializing", ncpu, "chains..."))

    fit = stan( 
        file = stanfile
        , fit = NA
        , data = stanInput
        , iter = n_iter
        , warmup = n_warmup
        , chains = ncpu
        , cores = ncpu
        , seed = runif(1)*10^6
        , refresh = 100
        , thin = thin_factor
        , init = lapply(1:ncpu, function(x){
            list(
                lambda = runif(length(Times)-1, 0.03, 0.05)
            )
        })
    )
    saveRDS( fit, file = stanRDS)
} else {
    message(paste("Fit exists:",stanRDS))
    message("Loading existing fit...")
    fit = readRDS(stanRDS)
}




    #   Summarize/generate quantities from the fits
    #   .............................

source("Scripts/foi/summarize/parameters.R")
source("Scripts/foi/summarize/susceptibility.R")
