#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Scatterplot of susceptibility status.")

parser$add_argument('outfile', help = 'Output file to save plot (relative to `outroot`).')
parser$add_argument('-outroot', default = '03-plots/susceptibility', help = 'Root directory to store plots.')

parser$add_argument('-obs.year.start', default = 1998, help = 'Year of 1st observation.')
parser$add_argument('-obs.year.end', default = 2019, help = 'Year of last observation.')

parser$add_argument('-S.x', required = T, help = 'Susceptibility states on x-axis.')
parser$add_argument('-S.y', required = T, nargs = "*", help = 'Susceptibility states on y-axis.')
parser$add_argument('-S.y.group', nargs = "*", help = 'Group (color) which `S.y` is associated with.')

parser$add_argument('-lab.x', default = 'Ground truth', help = 'Label of x-axis.')
parser$add_argument('-lab.y', default = 'Estimated', help = 'Label of y-axis.')



# parse arguments
inputArg = parser$parse_args()

outfile = with(inputArg, file.path(outroot, outfile))
dir.create(dirname(outfile), recursive = T)

S.breaks = c(0,1,2,Inf)
S.labels = c('Naive','Monotypic','Multitypic')
sumSusceptibility = function(x){
    x %>%
    mutate(S = cut(S
        , breaks = S.breaks
        , labels = S.labels
        , right = FALSE
    )) %>%
    group_by(year, age, S) %>%
    summarize(across(everything(), sum))
}

dat = 
    read_csv(inputArg$S.x) %>%
    sumSusceptibility %>%
    inner_join(
        mapply(function(f,g){
                read_csv(f) %>%
                sumSusceptibility %>%
                mutate(Source = g)
            }
            , f = inputArg$S.y
            , g = inputArg$S.y.group
            , SIMPLIFY = FALSE
        ) %>%
        do.call(what = rbind)
        , by = c('year', 'age', 'S')
    ) %>%
    mutate(Source = factor(Source, levels = unique(inputArg$S.y.group))) %>%
    ungroup %>%
    filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end))




    #   Compare susceptibility reconstructions
    #   ......................
    
compareValue = function(x,y){       
    coef = summary(lm(y ~ x + 0))$coefficients
    c(
        'Cor. coef.' = cor(x,y)
        , 'Ratio' = coef[1,1]
        , 'LB' = coef[1,1] - 1.96 * coef[1,2]
        , 'UB' = coef[1,1] + 1.96 * coef[1,2]
    ) %>%
    round(2) %>%
    formatC(digits = 2, format = "f", flag = "0")
}
plotCompareValueDf = function(textsize = 8){
    mat =
        split(dat, dat$Source) %>%
        lapply(function(dSource){
            out = rbind( ''
                , split(dSource, dSource$S) %>%
                    lapply(function(x){ with(x, compareValue(Val.x, Val.y))}) %>%
                    do.call(what = cbind)
            )
            rownames(out)[1] = as.character(dSource$Source[1])
            out    
        }) %>%
        do.call(what = rbind)

    ggtexttable( mat
        , theme = ttheme("blank", base_size = textsize, padding = unit(c(0.8, 0.8), "mm"))
    )
}






    #   Plot susceptibility scatterplot
    #   ...................
    
source("Scripts/configs/Aesthetics.R")

gScatter = 
    dat %>%
    ggplot(aes(x = Val.x, y = Val.y, color = Source))+
    geom_abline(slope = 1, linetype = 3, linewidth = 0.3)+
    geom_point(size = 1, shape = 1, stroke = 0.2)+
    facet_grid(. ~ S)+
    coord_fixed(ratio = 1)+
    # scale_color_brewer(palette = "Dark2")+
    scale_color_manual(values = Colors$Source)+
    xlab(inputArg$lab.x)+
    ylab(inputArg$lab.y)+
    xlim(c(0,1))+
    ylim(c(0,1))+
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )


gCorr = plotCompareValueDf(10)

g =
    ggarrange(
        gScatter %+% theme(legend.position = 'bottom')
        , gCorr
        , ncol = 2
        , nrow = 1
    )

ggsave(g
    , filename = outfile
    , width = 8
    , height = 2.5
)


# prop naive at age 9 in 2019
message("Proportion naive at age 9 in 2019 ....")
dat %>%
    filter(year == 2019, age == 9, S == "Naive") %>%
    group_by(Source) %>%
    summarize_at( vars(matches('[vV]al\\.[xy]$'))
        , .funs = ~ round(mean(.),2)
    )
