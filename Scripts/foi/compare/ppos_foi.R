#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Compare FOI from various P(pos) settings against a reference FOI series.")

parser$add_argument('outfile', help = 'Output file to save plot (relative to `outroot`).')
parser$add_argument('-outroot', default = '03-plots/ppos_foi', help = 'Root directory to store plots.')

parser$add_argument('-obs.year.start', default = 1998, help = 'Year of 1st observation.')
parser$add_argument('-fitdir.seroinc', default = '', help = 'Root directory of seroincidence fits with various P(pos) settings.')
parser$add_argument('-fitdir.seroprev', default = '', help = 'Root directory of seroprevalence fits with various P(pos) settings.')
parser$add_argument('-fitdir.ref', required = T, help = 'Directory with reference FOI series.')

# ppos from joint serology model
parser$add_argument('-ppos0', type = 'numeric')
parser$add_argument('-ppos1.long', type = 'numeric')


# parse arguments
inputArg = parser$parse_args()



outfile = with(inputArg, file.path(outroot, outfile))
dir.create(dirname(outfile), recursive = T)

source("Scripts/configs/Aesthetics.R")



    #   Import all FOI series
    #   .....................

# reference FOI series
dat.ref = 
    file.path(inputArg$fitdir.ref, "foi.csv") %>%
    read_csv %>%
    filter(year >= inputArg$obs.year.start) %>%
    select(year, val)
  
# function to characterize congruence with reference series
getCongruence = function(fitdir){
    list.dirs(fitdir, full.names = F, recursive = F) %>%
    lapply(function(x){
        ppos = str_split(x, '_') %>% unlist %>% as.numeric
        
        dat = 
            file.path(fitdir, x, "foi.csv") %>%
            read_csv(col_types = "dddi") %>%
            select(year, val) %>%
            inner_join(
                dat.ref
                , by = c('year')
                , suffix = c('', '.ref')
            )
            
        tibble(
            ppos0 = ppos[1]
            , ppos.long = ppos[2]
            , ppos.short = ppos[3]
            , Cor = with(dat, cor(val, val.ref))
            , Ratio = coefficients(lm(val ~ val.ref + 0, data = dat))
        )  
    }) %>%
    do.call(what = rbind)
}

# function to plot the congruence summaries
plotCongruence = function(dat, fillvar, midVal = 1, Trans = 'identity'){
    dat %>%
        ggplot(aes(x = ppos0, y = ppos.long))+
        geom_tile(aes(fill = {{fillvar}} ))+
        
        # mark P(pos) inferred from joint serology model
        geom_point(
            x = inputArg$ppos0
            , y = inputArg$ppos1.long
            , size = 3
            , stroke = 2
            , shape = 4
        )+
        
        scale_x_continuous("P(+|S=0)"
            , expand = c(0,0)
            , breaks = function(x) seq(0, 0.1, by = 0.02)
        )+
        scale_y_continuous("Long-term P(+|S=1)"
            , expand = c(0,0)
            , breaks = function(x) seq(0.9, 1, by = 0.02)
        )+
        facet_grid(Source ~ .)+
        scale_fill_gradient2(
            trans = Trans
            , low = "#272342"
            , mid = "#eeeeee"
            , high = "red"
            , midpoint = midVal
        )+
        theme(
            legend.position = "top"
            , axis.text.x = element_text(angle = 40, hjust = 1)
        )
}




dat = 
    rbind(
        getCongruence(inputArg$fitdir.seroprev) %>% mutate(Source = "Seroprevalence")
        , getCongruence(inputArg$fitdir.seroinc) %>% mutate(Source = "Seroincidence")
    )

g =
    ggarrange(
        plotCongruence(dat, Cor, midVal = 0.5)
        , plotCongruence(dat, Ratio, midVal = 0, Trans = 'log2')
        , nrow = 1
        , ncol = 2
        , align = 'hv'    
    )
    
ggsave(g
    , filename = outfile
    , width = 4
    , height = 3.5
)