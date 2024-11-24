 #!/usr/bin/env Rscript --vanilla

library(argparse)
library(ggpubr)
library(tidyverse); theme_set(theme_classic())
library(rstan)

parser = ArgumentParser(description = "Fit FOI from seroconversion data.")


parser$add_argument('-outbase', required = T)
parser$add_argument('-stanver', default = "time")

# model specification arguments
parser$add_argument('-Kappa', action = "store_true", default = FALSE)
parser$add_argument('-cohortKappa', action = "store_true", default = FALSE)
parser$add_argument('-cohortSigma', action = "store_true", default = FALSE)
parser$add_argument('-kappaEpochBreak', type = "integer")

# serology related arguments
parser$add_argument('-indir.sero', required = T)
parser$add_argument('-titer.threshold', nargs = "*", type = "double", default = c(10,20))
parser$add_argument('-ageMin.sero', type = "double", default = 3)
parser$add_argument('-ageMax.sero', type = "double", default = 30)

# infection related arguments (affects both serology and case)
parser$add_argument('-break.age', default = "Scripts/foi/configs/ageBreaks_cohort.txt")
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


stanfile = file.path("Scripts/foi/jointsero", stanver,"model.stan")

iFold = Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.integer
outdir = file.path(outbase, stanver, iFold)
dir.create(outdir, recursive = T)
stanRDS = file.path(outdir, "fit.RDS")

# load necessary functions
source("Scripts/configs/encoding_functions.R")
source("Scripts/configs/titer_functions.R")


getGammaShapeRate = function(m,v){
    rate = m/v
    shape = rate * m
    list('shape' = shape, 'rate' = rate, 'mean' = m, 'var' = v)
}
# priors for antibody kinetics parameters
Gamma = getGammaShapeRate(Gamma, Gamma_var)



    #   Data
    #   ....

source("Scripts/foi/prepare/jointsero.R")

ageBreaks = readLines(break.age)
ageBreaks = c(0, as.numeric(ageBreaks), Inf) %>% unique

yearBreaks = readLines(break.year)
yearBreaks = c(as.numeric(yearBreaks), Inf)

if(is.null(kappaEpochBreak)){ kappaEpochBreak = integer(0) }
kappaEpochBreakMat = cbind(
    1
    , c(which(Times %in% kappaEpochBreak) )
    , c(length(Times) - 1                 )
)




    #   Fit model
    #   .........
Studies = c("kps1", "kps2", "kps3", "kfcs")

dat.prev =
    dat.prev %>%
    mutate(study = factor(study, levels = Studies)) %>%
    arrange(study)
dat.inc =
    dat.inc %>%
    mutate(study = factor(study, levels = Studies)) %>%
    arrange(study)


if(!file.exists(stanRDS) | overwriteRDS){

    # ... seroincidence
    # presence matrices, for each individual (interval x time)
    # encoding time spent in the time interval
    Presence = 
            # for each individual
            mapply(function(timeBirth, timeStart, timeEnd){
                # Kappa %*% #### commented out because this model does not give away Kappa(a)
                            #### but estimates it from the data instead
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
                    , ageStart = head(ageBreaks, -1)
                    , ageEnd = tail(ageBreaks, -1)
                    , SIMPLIFY = FALSE
                ) %>%
                do.call(what = rbind))
            }
            , timeBirth = with(dat.inc, year1 - age1)
            , timeStart = with(dat.inc, year1)
            , timeEnd = with(dat.inc, year2)
            , SIMPLIFY = FALSE
        )
    
    prePresence = 
            # for each individual
            mapply(function(timeBirth, timeStart, timeEnd){
                # Kappa %*% #### commented out because this model does not give away Kappa(a)
                            #### but estimates it from the data instead
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
                    , ageStart = head(ageBreaks, -1)
                    , ageEnd = tail(ageBreaks, -1)
                    , SIMPLIFY = FALSE
                ) %>%
                do.call(what = rbind))
            }
            , timeBirth = with(dat.inc, year0)
            , timeStart = with(dat.inc, year0)
            , timeEnd = with(dat.inc, year1)
            , SIMPLIFY = FALSE
        )


    # ... seroprevalence
    # presence matrices (age bin x time) for each individual
    # encoding time spent in the interval
    Presence_prev =
        # for each individual
        mapply(function(timeBirth, timeObs){
                # Kappa %*% #### commented out because this model does not give away Kappa(a)
                            #### but estimates it from the data instead
                # calculate average Kappa of each time interval
                (mapply(function(ageStart, ageEnd){
                        out = encodeTimeInterval(
                            timeBirth + ageStart
                            , min(timeBirth + ageEnd, timeObs)
                            , Times
                        )
                        out[out < 0] = 0
                        out
                    }
                    , ageStart = head(ageBreaks, -1)
                    , ageEnd = tail(ageBreaks, -1)
                    , SIMPLIFY = FALSE
                ) %>%
                do.call(what = rbind))
            }
            , timeBirth = with(dat.prev, year - age)
            , timeObs = with(dat.prev, year)
            , SIMPLIFY = FALSE
        )


    # form input for stan
    stanInput = list(
        num_epoch = length(Times) - 1
        , num_epoch_kappa = ncol(kappaEpochBreakMat) - 1
        , break_epoch_kappa = kappaEpochBreakMat
        , num_thresh = length(titer.threshold)
        , num_study = Studies %>% length
        , cohortKappa = as.integer(cohortKappa)
        , cohortSigma = as.integer(cohortSigma)
        # seroincidence
        , num_intervals = nrow(dat.inc)
        , num_intervals_study = factor(dat.inc$study, levels = Studies) %>% table
        , Presence = Presence
        , prePresence = prePresence
        , obs_inc = dat.inc %>% 
            select(starts_with("pos2")) %>%
            lapply(as.integer) %>%
            do.call(what = cbind)
        , include_interval = dat.inc %>% 
            select(starts_with("pos1")) %>%
            lapply(function(x) as.integer(!x)) %>%
            do.call(what = cbind)
        # seroprevalence
        , num_indiv = nrow(dat.prev)
        , num_indiv_study = factor(dat.prev$study, levels = Studies) %>% table
        , Presence_prev = Presence_prev
        , obs_prev = dat.prev %>% 
            select(starts_with("pos")) %>%
            lapply(as.integer) %>%
            do.call(what = cbind)
        
        , Omega_long_shape_prior = Omega_shape
        , Omega_long_rate_prior = Omega_rate
        , Omega_short_shape_prior = Gamma$shape
        , Omega_short_rate_prior = Gamma$rate

        , Sigma_prior = Sigma
        , titer_thresh = adjustTiter(titer.threshold) %>% t
    )


    # Kappa(a)
    if(Kappa){
        stanInput$len_binned_kappa = length(ageBreaks) - 1
    } else {
        stanInput$len_binned_kappa = 1
    }

    
    set.seed(ncpu*190199)
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
                , binned_kappa_rel = matrix(1, stanInput$num_epoch_kappa, stanInput$len_binned_kappa - 1)
                , Omega_long = Omega_shape/Omega_rate
                , Omega_short_ref = Gamma$shape/Gamma$rate
                , Omega_short_mult = rep(1, stanInput$num_study-1)
                , Omega0_rel = 0.01
                , Sigma_ref = Sigma * runif(1, 0.9, 1.1)
                , psi = 100
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

# plot trace of log-posterior
ggsave(traceplot(fit, "lp__") + ggtitle(outdir)
    , filename = file.path(outdir, "lp.pdf")
    , width = 5
    , height = 4
)

source("Scripts/foi/summarize/parameters.R")
source("Scripts/foi/summarize/susceptibility_kappaEst.R")




    #   Goodness-of-fit
    #   ...............

source("Scripts/configs/Aesthetics.R")

waldInterval = function(x, n, conf.level = 0.95){
    p = x/n
    sd = sqrt(p*((1-p)/n))
    Z = qnorm(c(
        (1-conf.level) / 2
        , 1 - (1-conf.level)/2
    ))
    names(Z) = c('lowerVal', 'upperVal')
    lapply(Z, function(z){
        p + z*sd
    }) %>%
    do.call(what = cbind) %>%
    as_tibble
}

ageCuts = c(0, seq(8,11), 16, 20, 30, Inf)
cutAge = function(x, Cuts = ageCuts){
    cut(x, breaks = Cuts, include.lowest = T)
}

# extract posterior of expected values
est = extract(fit, paste0("expect_", c('inc','prev')))

# compute seroprevalence by age group, study
compare_prev = 
    dat.prev %>%
    ungroup %>%
    mutate(Age = cutAge(age)) %>%
    select(Age, study) %>%
    mutate(iCol = seq_along(Age)) %>%
    nest(data = iCol) %>%
    mutate(data = lapply(data, function(x){
        est$expect_prev[ ,x$iCol, ,drop = F] %>%
        apply(3, function(d){
            out =
                quantile(
                    rowSums(d) / ncol(d)
                    , c(0.025, 0.5, 0.975)
                ) %>%
                t %>%
                as_tibble
            colnames(out) = c('lowerVal', 'Val', 'upperVal')
            out
        }) %>%
        do.call(what = rbind) %>%
        mutate(Threshold = titer.threshold)
    })) %>%
    unnest(cols = data) %>%
    mutate(Type = "Expected") %>%
    rbind(
        dat.prev %>%
        mutate(Age = cutAge(age)) %>%
        select(Age, study, matches('^pos')) %>%
        pivot_longer(matches('^pos'), names_to = "Threshold", values_to = "pos") %>%
        group_by(Age, study, Threshold) %>%
        summarize(
            num_pos = sum(pos)
            , num_total = n()
        ) %>%
        filter(num_total > 5) %>%
        mutate(
            Val = num_pos/num_total
            , Wald = mapply(waldInterval, x = num_pos, n = num_total, SIMPLIFY = F)
        ) %>%
        unnest(cols = Wald) %>%
        select(-num_pos, -num_total) %>%
        mutate(
            Threshold = titer.threshold[as.integer(factor(Threshold, levels = unique(Threshold)))]
            , Type = "Observed"
        )
    )

# compute seroincidence by age group, study
compare_inc =
    dat.inc %>%
    ungroup %>%
    mutate(Age = cutAge(age1)) %>%
    select(Age, study, starts_with('pos1')) %>%
    mutate(iCol = seq_along(Age)) %>%
    nest(data = -c(Age, study)) %>%
    mutate(data = lapply(data, function(x){
        seq_along(titer.threshold) %>%
        lapply(function(iThresh){
            out =
                quantile(
                    rowSums(est$expect_inc[ ,x$iCol[!x[[iThresh]]], iThresh]) / sum(!x[[iThresh]])
                    , c(0.025, 0.5, 0.975)
                    , na.rm = T
                ) %>%
                t %>%
                as_tibble
            colnames(out) = c('lowerVal', 'Val', 'upperVal')
            out            
        }) %>%
        do.call(what = rbind) %>%
        mutate(Threshold = titer.threshold, Type = "Expected")
    })) %>%    
    unnest(cols = data) %>%
    rbind(
        dat.inc %>%
        mutate(Age = cutAge(age1)) %>%
        group_by(Age, study) %>%
        select(Age, study, starts_with('pos')) %>%
        nest(data = -c(Age, study)) %>%
        mutate(data = lapply(data, function(x){
            seq_along(titer.threshold) %>%
            lapply(function(iThresh){
                x %>%
                select(matches(paste0("\\.",iThresh,"$"))) %>%
                filter_at(vars(starts_with('pos1')), function(x) !x) %>%
                summarize_at(vars(starts_with('pos2')), list(num_pos = sum, num_total = length))
            }) %>%
            do.call(what = rbind) %>%
            mutate(Threshold = titer.threshold, Type = "Observed")
        })) %>%
        unnest(cols = data) %>%
        filter(num_total > 5) %>%
        mutate(
            Val = num_pos/num_total
            , Wald = mapply(waldInterval, x = num_pos, n = num_total, SIMPLIFY = F)
        ) %>%
        unnest(cols = Wald) %>%
        select(-num_pos, -num_total)
    )


plotCompare = function(x, yLabel){
    g =
        x %>%
        mutate(Threshold = factor(Threshold), Type = factor(Type)) %>%
        ggplot(aes(x = as.numeric(Age) + as.numeric(Threshold)/4 + as.numeric(Type)/8 - 0.5, y = Val, color = Threshold))+
        geom_linerange(aes(ymin = lowerVal, ymax = upperVal))+
        geom_point(aes(shape = Type, size = Type), stroke = 1.2)+
        scale_color_manual(values = Colors$Titer)+
        scale_shape_manual(values = c(16, 2))+
        scale_size_manual(values = c(1,2))+
        scale_x_continuous( "Age"
            , breaks = seq_along(levels(x$Age))
            , labels = levels(x$Age)
        )+
        ylab(yLabel)+
        facet_grid(. ~ study)+
        theme(
            axis.text.x = element_text(angle = 40, hjust = 1)
            , legend.spacing = unit(0, 'lines')
            , legend.key.height = unit(0.8, 'lines')
        )
    return(g)
}

g =
    ggarrange(
        plotCompare(compare_prev, "Seroprevalence")
        , plotCompare(compare_inc, "Seroincidence")
        # ,     gCase
        , nrow = 2
        , ncol = 1
        , heights = c(1,1)  
        , labels = 'auto'  
    )

ggsave(g
    , filename = file.path(outdir, "expect.pdf")
    , width = 6.5
    , height = 8
)



    #   Parameter estimates
    #   ...................

getSummary = function(Pattern){
    file.path(outdir,"est.csv") %>%
    read_csv %>%
    filter(grepl(Pattern, param))
    
}


    # ... Titer-related parameters ...

x =
    getSummary("^p_pos\\[") %>%
    mutate(iParam = 
        str_extract(param, '[0-9,]+') %>%
        strsplit(',') %>%
        lapply(function(x){
            x = as.integer(x)
            tibble(
                Study = Studies[x[1]]
                , Term = c('Long-term', 'Short-term')[x[2]]
                , iSuscept = c('S=0','S=1','S>1')[x[3]]
                , Threshold = titer.threshold[x[4]]
            )
        })
    ) %>%
    unnest(cols = iParam) %>%
    mutate(Study = factor(Study, levels = Studies)) %>%
    filter(iSuscept != 'S>1')

x %>%
    ggplot(aes(x = Study, color = factor(Threshold)))+
    geom_pointrange(aes(y = val, ymin = lowerVal, ymax = upperVal), shape = 1)+
    facet_grid(iSuscept ~ Threshold + Term, scale = 'free_y')+
    ylab('Prob. test positive')+
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1)
    )

x %>%
    select(Term, iSuscept, Threshold, val, lowerVal, upperVal) %>% 
    unique %>% 
    as.data.frame
    
    
# Average titer at measurement
getSummary("^Omega_long") %>%
    mutate(Study = 'All', Term = 'Long-term') %>%
    rbind(
        getSummary("^Omega_short\\[") %>%
        mutate(Study = Studies, Term = 'Short-term')        
    ) %>%
    mutate(Study = factor(Study, levels = c(Studies, 'All'))) %>%
    ggplot(aes(x = Study, color = Term))+
    geom_hline(yintercept = 0, linetype = 3)+
    geom_pointrange(aes(y = val, ymin = lowerVal, ymax = upperVal), shape = 1)+
    ylab('Adjusted titer')


# Assay SD
getSummary("^Sigma\\[") %>%
    mutate(Study = factor(Studies, levels = Studies)) %>%
    ggplot(aes(x = Study))+
    geom_pointrange(aes(y = val, ymin = lowerVal, ymax = upperVal), shape = 1)+
    ylab('Assay SD')


# CXR titer in naives
getSummary("^Omega0_rel") %>%
    mutate(Study = factor(Studies, levels = Studies)) %>%
    ggplot(aes(x = Study))+
    geom_pointrange(aes(y = val, ymin = lowerVal, ymax = upperVal), shape = 1)+
    ylab('CXR titer (relative to long-lived DENV monotypic titer)')

