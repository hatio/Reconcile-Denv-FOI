#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(rstan)
library(abind)

parser = ArgumentParser(description = "Combine posteriors across folds.")

parser$add_argument('indir')


# parse arguments
inputArg = parser$parse_args()



    #   Import posteriors
    #   .................
    
Params = c(
    # infection risk
    'lambda'
    , 'binned_kappa'
    # observation process: case
    , 'binned_phi_a_rel'
    , 'phi_t'
    , 'p_dhf_rel'
    # observation process: serology
    , 'Omega_short'
    , 'Omega_short_mult'
    , 'Omega_long'
    , 'Omega0_rel'
    , 'Sigma'
    , 'p_pos'
)

est = 
    inputArg$indir %>%
    list.dirs(recursive = F, full.names = F) %>%
    lapply(function(iFold){
        fit =
            file.path(inputArg$indir, iFold, "fit.RDS") %>% 
            readRDS %>%
            extract
        fit[intersect(names(fit), Params)]
    })



    #   Export combined summaries of estimates
    #   .........................

summarizeParam = 
    function(param){

        x =
            lapply(est, function(x) x[[param]]) %>%
            abind(along = 1)

        # if parameter is not estimated
        if(is.null(x)){ return(NULL) }
        # if parameter is a single number
        if(length(dim(x)) == 1){
            x = as_tibble_row(quantile(x, c(0.025, 0.5, 0.975)))
            colnames(x) = c('lowerVal', 'val', 'upperVal')
            return(x %>% mutate(param))
        }
        # else
        x =    
            apply(x
                , seq_along(dim(x))[-1]
                , quantile
                , c(0.025, 0.5, 0.975)
            ) %>%
            as_tibble %>%
            t %>%
            as.data.frame
            
        colnames(x) = c('lowerVal', 'val', 'upperVal')
        x %>%
            mutate(param = paste0(param, '[', gsub('\\.', ',', gsub('^V', '', rownames(x))), ']')) %>%
            as_tibble
    }



#   FOI
#   ...

Years =
    inputArg$indir %>%
    list.dirs(recursive = F, full.names = F) %>%
    lapply(function(iFold){
        file.path(inputArg$indir, iFold, "foi.csv") %>% 
        read_csv %>%
        select(year)
    }) %>%
    unique

summarizeParam("lambda") %>%
    select(-param) %>%
    cbind(Years[[1]]) %>%
    write_csv(file.path(inputArg$indir, 'foi.csv'))


#   Other parameters
#   ................

fitsum = 
    Params %>%
    setdiff("lambda") %>%
    lapply(summarizeParam) %>%
    do.call(what = rbind)

if(!is.null(fitsum)){
    fitsum %>%
    write_csv(file.path(inputArg$indir, 'est.csv'))    
}
