#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())

parser = ArgumentParser(description = "Format estimated FOI for publication.")

parser$add_argument('indir', help = 'Directory of the FOI estimates.')


# parse arguments
inputArg = parser$parse_args()



source("Scripts/configs/general_functions.R")

    #   Import
    #   ......
    
file.path(inputArg, "foi.csv") %>% 
    read_csv %>%
    mutate_at(vars(matches('[vV]al$')), formatValue) %>%
    mutate(
        FOI = paste0(val, " (", lowerVal, ", ", upperVal, ")")
    ) %>%
    select(year, FOI) %>%
    write_csv(
        file.path(inputArg$indir, "foi_formatted.csv")
    )
    