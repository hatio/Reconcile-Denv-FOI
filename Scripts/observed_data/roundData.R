
library(tidyverse)

indir = "01-processedData/serology/gmt_orig"
outdir = "01-processedData/serology/gmt"
dir.create(outdir, recursive = T)

list.files(indir) %>%
    lapply(function(f){
        file.path(indir, f) %>%
        read_csv %>%
        mutate(
            dateCollect = 
                paste(
                    format(dateCollect, "%Y-%m")
                    # , ceiling(as.integer(format(dateCollect, "%m")) / 3) * 3
                    , "01"
                    , sep = "-"
                ) %>% as.Date
            , age = round(age)
        ) %>%
        nest(-id) %>%
        mutate(id = seq_along(id)) %>%
        unnest(cols = data) %>%
        write_csv(file.path(outdir, f))
    }) %>%
    invisible

