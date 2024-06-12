#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Plot distribution of infection risk by age (Kappa).")

parser$add_argument('outfile', help = 'Output file to save plot (relative to `outroot`).')
parser$add_argument('-outroot', default = '03-plots/kappa', help = 'Root directory to store plots.')

parser$add_argument('-truthdir', default = '', help = "Directory storing true Tau(t), Kappa(a).")
parser$add_argument('-fitdir.seroinc', default = '', help = 'Directory of fits from seroincidence data.')
parser$add_argument('-fitdir.seroprev', default = '', help = 'Directory of fits from seroprev data.')
parser$add_argument('-fitdir.case', default = '', help = 'Directory of fits from case count data.')
parser$add_argument('-fitdir.joint', default = '', help = 'Directory of joint fits.')
parser$add_argument('-fitdir.jointsero', default = '', help = 'Directory of joint serology fits.')
parser$add_argument('-break.age.sero', default = "Scripts/foi/configs/ageBreaks_cohort.txt")
parser$add_argument('-break.age.case', default = "Scripts/foi/configs/ageBreaks.txt")

parser$add_argument('-ref', nargs = "*", default = 'Truth', help = 'Which Kappa series to treat as reference(s): no CI, fainted color.')
parser$add_argument('-splitFold', action="store_true", default=FALSE, help = 'Whether to show estimates of each fold separately.')



# parse arguments
inputArg = parser$parse_args()

outfile = with(inputArg, file.path(outroot, outfile))
dir.create(dirname(outfile), recursive = T)


source("Scripts/configs/Aesthetics.R")
source("Scripts/configs/general_functions.R")

# age breaks of Kappa fits
ageBreaks.sero = readLines(inputArg$break.age.sero) %>% as.integer
ageBreaks.sero = c(ageBreaks.sero, max(ageBreaks.sero)+3)

ageBreaks.case = readLines(inputArg$break.age.case) %>% as.integer
ageBreaks.case = c(ageBreaks.case, max(ageBreaks.case)+3)


# function to export inferred age-specific risk for later use
exportKappa = function(f, outstem){
    f %>%
        mutate(
            ageEnd = ifelse(is.na(ageEnd), 100, ageEnd)
            , Age = mapply(seq, ageStart, ageEnd-1, SIMPLIFY = F)
            , Kappa = val
        ) %>%
        select(Age, Kappa) %>%
        unnest(cols = Age) %>%
        write_csv(file = file.path(outstem, "Kappa.csv"))   
}
exportKappaFormatted = function(f, outstem){
    f %>%
        mutate_at(vars(matches('[vV]al$')), formatValue) %>%
        mutate(
            Age = paste0(ageStart, "-", ageEnd-1, " yrs")
            , Kappa = paste0(val, " (", lowerVal, ", ", upperVal, ")")
        ) %>%
        select(Age, Kappa) %>%
        write_csv(file = file.path(outstem, "Kappa_formatted.csv"))   
}
# function to process kappa summaries from summaries of estimates
processKappaSummary = function(f, ageBreaks){
    f =
        f %>%
        filter(grepl('^binned_kappa\\[', param))
    rbind(f
        , tail(f,1) %>% mutate_at(vars(matches('[vV]al$')), function(x) 1)
    ) %>%
    mutate(
        ageStart = ageBreaks[seq_along(param)]
        , ageEnd = ageBreaks[seq_along(param)+1]
    ) %>%
    select(-param)
}


# kappa: from seroprevalence
if(inputArg$fitdir.seroprev != ''){
    if(inputArg$splitFold){
        dat.prev = 
            inputArg$fitdir.seroprev %>%
            list.dirs(full.names = F, recursive = F) %>%
            Filter(f = function(x){ grepl('^[0-9]+$', x) }) %>%
            lapply(function(iFold){
                f = file.path(inputArg$fitdir.seroprev, iFold, "est.csv")
                if(!file.exists(f)){ return(NULL) }
                f = f %>%
                    read_csv(col_types = "cddddd") %>%
                    processKappaSummary(ageBreaks = ageBreaks.sero) %>%
                    mutate(
                        iFold = iFold
                        , .before = 1
                    )
                # export inferred age-specific risk for later use
                exportKappa(f, file.path(inputArg$fitdir.seroprev, iFold))
                return(f)                
            }) %>%
            do.call(what = rbind)
    } else {
        dat.prev =
            file.path(inputArg$fitdir.seroprev, "est.csv") %>%
            read_csv(col_types = "dddc") %>%
            processKappaSummary(ageBreaks = ageBreaks.sero) %>%
            mutate(iFold = "combined")
        exportKappa(dat.prev, file.path(inputArg$fitdir.seroprev))
    }
}


# kappa: from seroincidence
if(inputArg$fitdir.seroinc != ''){
    dat.inc =
        file.path(inputArg$fitdir.seroinc, 'est.csv') %>%
        lapply(function(f){
            if(!file.exists(f)){ return(NULL) }
            f = f %>%
                read_csv(col_types = "cddddd") %>%
                filter(grepl('^binned_kappa\\[', param)) 
            f = rbind(f
                , tail(f,1) %>% mutate_at(vars(matches('[vV]al$')), function(x) 1)
            ) %>%
                mutate(
                    ageStart = ageBreaks.sero[seq_along(param)]
                    , ageEnd = ageBreaks.sero[seq_along(param)+1]
                ) %>%
                select(-param) %>%
                mutate(
                    iFold = "all"
                    , .before = 1
                )
            # export inferred age-specific risk for later use
            exportKappa(f, inputArg$fitdir.seroinc)
            return(f)
        }) %>%
        do.call(what = rbind)
}

 
# kappa: from case data
if(inputArg$fitdir.case != ''){
    dat.case = 
        file.path(inputArg$fitdir.case, "est.csv") %>%
        read_csv(col_types = "cddddd") %>%
        filter(grepl('^binned_kappa_rel', param)) 
    dat.case =
        head(dat.case,1) %>% mutate_at(vars(matches('[vV]al$')), function(x) 1) %>%
        rbind(dat.case) %>%
        rbind(tail(dat.case,1) %>% mutate_at(vars(matches('[vV]al$')), function(x) 1)) %>%
        mutate(
            ageStart = ageBreaks.case[seq_along(param)]
            , ageEnd = ageBreaks.case[seq_along(param)+1]
        ) %>%
        select(-param) %>%
        mutate(iFold = "all", .before = 1)

    # export inferred age-specific risk for later use
    exportKappa(dat.case, inputArg$fitdir.case)
    # export formatted version for manuscript
    exportKappaFormatted(dat.case, inputArg$fitdir.case)
}


 


# kappa: Truth
if(inputArg$truthdir != ''){
    dat.truth = 
        file.path(inputArg$truthdir, 'Kappa.csv') %>% read_csv %>%
        filter(Age <= max(ageBreaks.case)) %>%
        mutate(
            iAgeGroup = c(T,head(Kappa,-1) != tail(Kappa,-1)) %>% cumsum
        ) %>%
        group_by(iAgeGroup) %>%
        summarize(
            ageStart = min(Age)
            , ageEnd = max(Age) + 1
            , val = unique(Kappa)
        ) %>%
        mutate(lowerVal = val, upperVal = val, iFold = "all")
}




# kappa: from joint model
if(inputArg$fitdir.joint != ""){
    
    if(inputArg$splitFold){
        dat.joint = 
            inputArg$fitdir.joint %>%
            list.dirs(full.names = F, recursive = F) %>%
            Filter(f = function(x){ grepl('^[0-9]+$', x) }) %>%
            lapply(function(iFold){
                f = file.path(inputArg$fitdir.joint, iFold, "est.csv")
                if(!file.exists(f)){ return(NULL) }
                f = f %>%
                    read_csv(col_types = "cddddd") %>%
                    processKappaSummary(ageBreaks = ageBreaks.case) %>%
                    mutate(
                        iFold = iFold
                        , .before = 1
                    )
                # export inferred age-specific risk for later use
                exportKappa(f, file.path(inputArg$fitdir.joint, iFold))
                return(f)  
            }) %>%
            do.call(what = rbind)
    } else {
        dat.joint =
            file.path(inputArg$fitdir.joint, "est.csv") %>%
            read_csv(col_types = "dddc") %>%
            processKappaSummary(ageBreaks = ageBreaks.case) %>%
            mutate(iFold = "combined")
        exportKappa(dat.joint, file.path(inputArg$fitdir.joint))
    }
}



# kappa: from jointsero model
if(inputArg$fitdir.jointsero != ""){
    
    if(inputArg$splitFold){
        dat.jointsero = 
            inputArg$fitdir.jointsero %>%
            list.dirs(full.names = F, recursive = F) %>%
            Filter(f = function(x){ grepl('^[0-9]+$', x) }) %>%
            lapply(function(iFold){
                f = file.path(inputArg$fitdir.jointsero, iFold, "est.csv")
                if(!file.exists(f)){ return(NULL) }
                f = f %>%
                    read_csv(col_types = "cddddd") %>%
                    processKappaSummary(ageBreaks = ageBreaks.sero) %>%
                    mutate(
                        iFold = iFold
                        , .before = 1
                    )
                # export inferred age-specific risk for later use
                exportKappa(f, file.path(inputArg$fitdir.jointsero, iFold))
                return(f)  
            }) %>%
            do.call(what = rbind)
    } else {
        dat.jointsero =
            file.path(inputArg$fitdir.jointsero, "est.csv") %>%
            read_csv(col_types = "dddc") %>%
            processKappaSummary(ageBreaks = ageBreaks.sero) %>%
            mutate(iFold = "combined")
        exportKappa(dat.jointsero, file.path(inputArg$fitdir.jointsero))
    }
    
    # export formatted version for manuscript
    exportKappaFormatted(dat.jointsero, inputArg$fitdir.jointsero)
}




Datasets = c(
        "dat.inc" = "Seroincidence"
        , "dat.prev" = "Seroprevalence"
        , "dat.case" = "Case"
        , "dat.joint" = "Joint"
        , "dat.jointsero" = "Joint sero."
        , "dat.truth" = "Truth"
    )
Datasets = Datasets[names(Datasets) %in% ls(pattern = '^dat\\.')]




# combine data sources
dat.all = 
    names(Datasets) %>%
    lapply(function(x){
        get(x) %>% 
        mutate(Source = Datasets[[x]]) %>%
        select(Source, iFold, ageStart, ageEnd, val, lowerVal, upperVal)
    }) %>%
    do.call(what = rbind) %>%
    mutate(Source = factor(Source, levels = Datasets))



# remove CI from reference series
dat.all =
    dat.all %>%
    mutate(
        isRef = Source %in% inputArg$ref
        , lowerVal = ifelse(isRef, NA, lowerVal)
        , upperVal = ifelse(isRef, NA, upperVal)
    )


g = 
    dat.all %>%
    ggplot()+
    geom_hline(yintercept = 1, linetype = 3, color = '#dddddd')+
    # bounds ....
    geom_rect(
        aes(xmin = ageStart, xmax = ageEnd, ymin = lowerVal, ymax = upperVal, fill = Source)
        , alpha = 0.15
    )+
    # averages ....
    geom_step(
        aes(
            x = ageStart, y = val
            , color = Source, group = paste(Source,iFold)
            , linewidth = isRef
        )
    )+
    scale_linewidth_manual(values = c(1, 0.4), guide = "none")+
    scale_color_manual(values = Colors$Source, aesthetics = c('color', 'fill'))+
    xlab("Age")+
    ylab("Infection risk\nrelative to youngest age group")+
    coord_cartesian(xlim = c(
        0
        , dat.all %>%
            filter(Source != "Truth") %>%
            with(max(ageEnd, na.rm = T))
    ))

ggsave(g
    , filename = outfile
    , width = 5
    , height = 3
)

