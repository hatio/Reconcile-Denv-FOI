#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(rstan)

parser = ArgumentParser(description = "Fit FOI from seroconversion data.")


parser$add_argument('-outbase', required = T)
parser$add_argument('-stanver', default = "time")

# case related arguments
parser$add_argument('-casedir', required = T)
parser$add_argument('-popfile', required = T)
parser$add_argument('-model.mod', default = 'm1h.rta')

# infection related arguments (affects both serology and case)
parser$add_argument('-break.age', default = "Scripts/foi/configs/ageBreaks.txt")
parser$add_argument('-break.year', default = "Scripts/foi/configs/yearBreaks.txt")
parser$add_argument('-T.delta', type = "double", default = 1)
parser$add_argument('-T.pre', type = "integer", default = 10, help = "Number of time-varying FOI pre-study period before assuming constant FOI.")

# MCMC configurations
parser$add_argument('-ncpu', type = "integer", default = 3)
parser$add_argument('-thin_factor', type = "integer", default = 5)
parser$add_argument('-n_iter', type = "integer", default = 2000)
parser$add_argument('-n_warmup', type = "integer", default = 200)
parser$add_argument('-overwriteRDS', action = "store_true", default = FALSE)



# parse arguments
inputArg = parser$parse_args()
attach(inputArg)


stanfile = file.path("Scripts/foi/casecount", stanver,"model.stan")


outdir = file.path(outbase, stanver, model.mod)
dir.create(outdir, recursive = T)
stanRDS = file.path(outdir, "fit.RDS")

# load necessary functions
source("Scripts/configs/encoding_functions.R")
source("Scripts/configs/general_functions.R")




    #   Data
    #   ....

source("Scripts/foi/prepare/casecount.R")
Times = Years

ageBreaks = readLines(break.age)
ageBreaks = c(0, as.numeric(ageBreaks), Inf) %>% unique

yearBreaks = readLines(break.year)
yearBreaks = c(as.numeric(yearBreaks), Inf) %>% unique


    #   Fit model
    #   .........

if(!file.exists(stanRDS) | overwriteRDS){

    # form input for stan
    stanInput = list(
        num_epoch = length(Times) - 1
        # case data
        , len_cohorts = nrow(pop)
        , len_times_case = ncol(case)
        , len_Ag = nrow(ageGroups)
        , len_ages = ncol(ageGroups)
        , mapAge = ageGroups
        , obs_case = as.matrix(case)
        , pop = as.matrix(pop)
        
        , naVal = naVal
    )


    # ... case data
    # which epoch does each time index map to?
    stanInput$epoch_time = 
        cut( 1:stanInput$len_cohorts
            , breaks = c(
                0
                , seq(stanInput$len_cohorts - stanInput$len_times_case - T.pre*T.delta + 1, stanInput$len_cohorts, by = 1)
            ) %>%
            unique
        ) %>%
        as.integer
    stanInput$len_epoch = length(unique(stanInput$epoch_time))

    # Kappa(a)
    if(grepl("m3", model.mod)){
        stanInput$len_binned_kappa = length(ageBreaks) - 1
        stanInput$index_kappa = cut( 1:stanInput$len_ages, breaks = ageBreaks, right = FALSE) %>% as.integer
    } else {
        stanInput$len_binned_kappa = 1
        stanInput$index_kappa = rep(1, stanInput$len_ages)
    }

    # Phi(a)
    if(grepl("rt*a", model.mod)){
        stanInput$len_binned_phi_a = length(ageBreaks) - 1
        stanInput$index_phi_a = cut( 1:stanInput$len_ages, breaks = ageBreaks, right = FALSE) %>% as.integer
    } else {
        stanInput$len_binned_phi_a = 1
        stanInput$index_phi_a = rep(1, stanInput$len_ages)
    }

    # Phi(t)
    if(grepl("rt", model.mod)){
        stanInput$len_binned_phi_t = length(yearBreaks) - 1
        stanInput$index_phi_t = cut(max(Years) - (stanInput$len_times_case : 1), breaks = yearBreaks, right = FALSE) %>% as.integer
    } else {
        stanInput$len_binned_phi_t = 1
        stanInput$index_phi_t = rep(1, stanInput$len_times_case)
    }


    
    set.seed(ncpu*28945)
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
                # , binned_kappa_rel = rep(1, stanInput$len_binned_kappa - 1)
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
source("Scripts/foi/summarize/susceptibility_kappaEst.R")



    #   Format estimates for manuscript
    #   ...............................

Params = c(
    'p_dhf_rel'
    , 'binned_phi_t'
    , 'binned_phi_a_rel'
)
est = 
    file.path(outdir, "est.csv") %>%
    read_csv(col_types = "cddddd")

getFormattedParam = function(parname){
    est %>%
        filter(grepl(paste0('^', parname, '$|^', parname, '\\['), param)) %>%
        mutate_at(vars(matches('[vV]al$')), formatValue) %>%  
        mutate(value = paste0(val, " (", lowerVal, ", ", upperVal, ")")) %>%
        select(param, value)
}


# phi(t)
    
    tibble(
        yearStart = head(yearBreaks, -1)
        , yearEnd = tail(yearBreaks, -1) - 1
    ) %>%
    mutate(
        Year = paste0(yearStart, "-", yearEnd)
    ) %>%
    select(Year) %>%
    cbind(getFormattedParam('binned_phi_t')) %>%
    write_csv(file.path(outdir, "formatted_phi_t.csv"))
    
    
# phi(a)

    tibble(
        ageStart = head(ageBreaks, -1)
        , ageEnd = tail(ageBreaks, -1) - 1
    ) %>%
    # remove reference class
    tail(-1) %>%
    mutate(
        Age = paste0(ageStart, "-", ageEnd)
    ) %>%
    select(Age) %>%
    cbind(getFormattedParam('binned_phi_a_rel')) %>%
    write_csv(file.path(outdir, "formatted_phi_a.csv"))

# P(DHF)

    getFormattedParam('p_dhf') %>%
    write_csv(file.path(outdir, "formatted_p_dhf.csv"))



    #   Plot expected case counts
    #   .........................

# compute age distribution of cases by year
ageGroup.case =
    ageGroups %>%
    apply(1, function(x){
        which(x > 0) %>%
        unique %>%
        as.integer
    })
compare_case =
    extract(fit, c("expect_case"))$expect_case %>%
    apply(2:3, quantile, c(0.025, 0.5, 0.975)) %>%
    apply(1, function(x){ x

        1:ncol(case) %>%
        lapply(function(y){
            tibble(ageGroup = ageGroup.case, Val = x[ ,y]) %>%
            mutate(Year = colnames(case)[y])
        }) %>%
        do.call(what = rbind)
    }, simplify = F)
compare_case =
    mapply(function(x,newName){
            x %>%
            rename_with(function(a) newName, Val)
        }
        , x = compare_case
        , newName = c('lowerVal','Val','upperVal')
        , SIMPLIFY = F
    ) %>%
    Reduce(f = full_join) %>%
    mutate(Type = "Expected")

gCase = 
    compare_case %>%
    ggplot(aes(x = ageGroup, y = Val))+
    # observed
    geom_col(data = 
        colnames(case) %>%
            lapply(function(y){
                tibble(ageGroup = ageGroup.case, Val = case[[y]]) %>%
                filter(Val != -1) %>%
                mutate(Year = y)
            }) %>%
            do.call(what = rbind)
        , width = 1
    )+    
    # expected
    geom_ribbon(aes(ymin = lowerVal, ymax = upperVal), alpha = 0.4)+
    geom_line()+
    facet_wrap( ~ Year, nrow = 4)+
    ylab("Num. case")+
    xlab("Age")
    
ggsave(gCase
    , filename = file.path(outdir, "expect.pdf")
    , width = 6.5
    , height = 8
)
