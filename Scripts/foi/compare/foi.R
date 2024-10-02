#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse); theme_set(theme_classic())
library(cowplot)
library(ggpubr)

parser = ArgumentParser(description = "Plot FOI fits.")

parser$add_argument('outfile', help = 'Output file to save plot (relative to `outroot`).')
parser$add_argument('-outroot', default = '03-plots/foi', help = 'Root directory to store plots.')
parser$add_argument('-inset', action = "store_true", default = FALSE, help = 'Plot scatterplot as inset instead of side panel.')


parser$add_argument('-obs.year.start', default = 1998, help = 'Year of 1st observation.')
parser$add_argument('-obs.year.end', default = 2019, help = 'Year of last observation.')
parser$add_argument('-fitdir.case', default = '', help = 'Directory of fits from case count data.')
parser$add_argument('-fitdir.seroinc', default = '', help = 'Directory of fits from seroincidence data.')
parser$add_argument('-fitdir.seroprev', default = '', help = 'Directory of fits from seroprev data.')
parser$add_argument('-seroprev.allBleed', action="store_true", default=FALSE, help = 'Whether to use seroprevalence fits where all bleeds from every individual were used.')
parser$add_argument('-fitdir.joint', default = '', help = 'Directory of joint fits.')
parser$add_argument('-splitFold', action="store_true", default=FALSE, help = 'Whether to show estimates of each fold separately.')

parser$add_argument('-truthdir', default = '01-processedData/simcohort/m1/infection_param/mock', help = "Directory storing true Tau(t), Kappa(a) for simulations.")
parser$add_argument('-ref', default = 'Truth', help = 'Which FOI series to use as reference for comparisons.')
parser$add_argument('-dropRatioCi', action="store_true", default=FALSE, help = 'Drop CI of ratio compared to reference series (default: FALSE).')
parser$add_argument('-foiAxisMax', type = "double", help = 'Max of FOI axis.')
parser$add_argument('-foiAxisIncrement', type = "double", help = 'Increment of FOI axis.')

parser$add_argument('-salje', action="store_true", default=FALSE, help = 'Overlay FOI estimates from Salje, 2018.')



# parse arguments
inputArg = parser$parse_args()


# Compare incidence rates of 1st infections
# simple vs adjusted 1st infection incidence rates in simulated cohort
outfile = with(inputArg, file.path(outroot, outfile))
dir.create(dirname(outfile), recursive = T)

source("Scripts/configs/Aesthetics.R")



    #   Import all FOI series
    #   .....................

if(inputArg$ref == "Truth"){
    # foi: Simulation ground truth
    trueFoi = file.path(inputArg$truthdir, 'Tau.csv') %>% read_csv
    trueKappa = file.path(inputArg$truthdir, 'Kappa.csv') %>% read_csv
    trueFoi$Kappa = list(trueKappa %>% 
        filter(between(Age,5,15)) %>%
        mutate(
            AgeGroup = c(F, head(Kappa,-1)!=tail(Kappa,-1))
            , AgeGroup = cumsum(AgeGroup)
        ) %>%
        group_by(AgeGroup) %>%
        summarize(
            AgeGroup = paste('Age', min(Age), 'to', max(Age))
            , AgeGroup = factor(AgeGroup, levels = AgeGroup)
            , Kappa = unique(Kappa)
        )
    )
    dat.truth = trueFoi %>%
        # to make geom_step of last FOI year have width
        rbind(tail(trueFoi, 1) %>% mutate(Year = Year+1)) %>%
        rename(year = Year, val = Tau) %>%
        mutate(lowerVal = val, upperVal = val, iFold = "All") %>%
        filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end))
}


# foi from seroprevalence
if(inputArg$fitdir.seroprev != ""){

    if(inputArg$splitFold){
        dat.prev = 
            inputArg$fitdir.seroprev %>%
            list.dirs(full.names = F, recursive = F) %>%
            Filter(f = function(x){
                (x == "NA") == inputArg$seroprev.allBleed
            }) %>%
            lapply(function(iFold){
                f = file.path(inputArg$fitdir.seroprev, iFold, "foi.csv")
                if(!file.exists(f)){ return(NULL) }
                f %>%
                read_csv(col_types = "iddd") %>%
                mutate(
                    iFold = ifelse(is.na(iFold), "all", iFold)
                    , .before = 1
                )
            }) %>%
            do.call(what = rbind)
    } else {
        dat.prev =
            file.path(inputArg$fitdir.seroprev, "foi.csv") %>%
            read_csv(col_types = "dddi") %>%
            mutate(iFold = "combined")
    }
    
    dat.prev =
        dat.prev %>%
            mutate(
                iFold = factor(iFold)
            ) %>%
            filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end))
}


# foi: from seroincidence
if(inputArg$fitdir.seroinc != ""){
    dat.inc = file.path(inputArg$fitdir.seroinc, 'foi.csv') %>%
        lapply(function(f){
            if(!file.exists(f)){ return(NULL) }
            f %>%
                read_csv(col_types = "iddd") %>%
                mutate(
                    iFold = "all"
                    , .before = 1
                )
        }) %>%
        do.call(what = rbind) %>%
        filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end)) %>%
        filter(year %in% (
            file.path(inputArg$fitdir.seroinc, 'Times_data.txt') %>%
            readLines %>%
            as.numeric
        ))
}


# foi from case data
if(inputArg$fitdir.case != ""){

    dat.case = file.path(inputArg$fitdir.case, "foi.csv") %>%
        read_csv %>%
        filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end)) 
    dat.case = dat.case %>%
        # to make geom_ribbon of first FOI year match cohorts
        rbind(head(dat.case, 1) %>% mutate(year = inputArg$obs.year.start)) %>%
        unique %>%
        mutate(iFold = "All")
}


# joint foi
if(inputArg$fitdir.joint != ""){

    if(inputArg$splitFold){
        dat.joint =
            inputArg$fitdir.joint %>%
            list.dirs(full.names = F, recursive = F) %>%
            lapply(function(iFold){
                f = file.path(inputArg$fitdir.joint, iFold, "foi.csv")
                if(!file.exists(f)){ return(NULL) }
                f = 
                    f %>%
                    read_csv(col_types = "iddd") %>%
                    mutate(
                        iFold = ifelse(is.na(iFold), "all", iFold)
                        , .before = 1
                    )
            }) %>%
            do.call(what = rbind)         
    } else {
        dat.joint =
            file.path(inputArg$fitdir.joint, "foi.csv") %>%
            read_csv(col_types = "dddi") %>%
            mutate(iFold = "combined")
    }

    dat.joint =
        dat.joint %>%
        mutate(
            iFold = factor(iFold)
        ) %>%
        filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end))
}




    #   Plot FOI series
    #   ...............

plotSeries = function(g, xname){
    g +
    geom_rect(data = get(paste0('dat.', tolower(xname))) %>%
            group_by(iFold) %>%
            arrange(year) %>%
            mutate(yearEnd = c(tail(year,-1), tail(year,1)+1))
        , aes(xmin = year, xmax = yearEnd, ymin = lowerVal, ymax = upperVal, group = iFold)
        , alpha = 0.3
        , fill = Colors$Source[[xname]]
    )+
    geom_step(data = get(paste0('dat.', tolower(xname))) %>%
            group_by(iFold) %>%
            filter(year == max(year)) %>%
            mutate(year = year + 1) %>%
            rbind(get(paste0('dat.', tolower(xname)))) %>%
            arrange(year)
        , aes(y = val, group = iFold)
        , direction = 'hv'
        , color = Colors$Source[[xname]]
    )
}




# initiate plot with reference series
gTime = 
    ggplot(mapping = aes(x = year))+
    scale_color_grey()+
    ylab('Per-serotype FOI')+
    guides(fill = "none")+
    theme(
        axis.title.x = element_blank()
        , legend.position = c(1,1)
        , legend.justification = c(1,1)
    )
gTime = plotSeries(gTime, inputArg$ref)


# plot case (if not already plotted)
if(inputArg$ref != 'Case' & exists('dat.case')){
    gTime = plotSeries(gTime, 'Case')
}

# plot joint (if not already plotted)
if(inputArg$ref != 'Joint' & exists('dat.joint')){
    gTime = plotSeries(gTime, 'Joint')
}


# plot: foi from seroprevalence
if(exists('dat.prev')){
    gTime = 
        gTime +
        geom_rect(data = dat.prev %>%
                group_by(iFold) %>%
                mutate(yearEnd = c(tail(year,-1), tail(year,1)+1))
            , aes(xmin = year, xmax = yearEnd, ymin = lowerVal, ymax = upperVal)
            , fill = "#F7DD5C", alpha = 0.3
        )+
        geom_boxplot(data = dat.prev %>% mutate(year = year + 0.5)
            , aes(y = val, group = year)
            , color = Colors$Source[['Seroprevalence']]
        )
}

# plot: foi from seroincidence
if(exists('dat.inc')){
    gTime = 
        gTime +
        geom_linerange(data = dat.inc %>% mutate(year = year + 0.5)
            , aes(ymin = lowerVal, ymax = upperVal), color = '#EB4134')+
        geom_point(data = dat.inc %>% mutate(year = year + 0.5)
            , aes(y = val), size = 0.7, color = '#EB4134')

}


# plot: overlay estimates from Salje, 2018
if(inputArg$salje){
    gTime = 
        gTime +
        geom_pointrange(data = 
                read_csv("01-processedData/literature/salje.csv") %>%
                filter(between(year, inputArg$obs.year.start, inputArg$obs.year.end)) %>%
                mutate(year = year + 0.5)
            , aes(y = foi, ymin = foi_lower, ymax = foi_upper)
            , shape =1, size = 0.3, stroke = 0.5
        )
    
}






    #   Compare FOI series
    #   ..................


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
compareValueDf = function(x){
    x %>%
    inner_join(
        get(paste0('dat.', tolower(inputArg$ref))) %>%
            rename_at(vars(matches('[vV]al$')), function(x) paste0(x,'.ref'))
        , by = 'year'
    ) %>%
    with(compareValue(val.ref, val))    
}
plotCompareValueDf = function(textsize = 8){
    Sources =
        ls(pattern = '^dat\\.', name = .GlobalEnv) %>%
        gsub(pattern = '^dat\\.', replacement = '') %>%
        setdiff('truth') %>%
        setdiff(tolower(inputArg$ref)) 
    print(Sources)
    mat =
        Sources %>%
        lapply(function(x) get(paste0('dat.', x))) %>%
        lapply(compareValueDf) %>%
        do.call(what = cbind)
    mat = rbind(c(
            'prev' = 'Seroprev.'
            , 'inc' = 'Seroinc.'
            , 'case' = 'Case'
            , 'joint' = 'Joint'
        )[Sources], mat)
    colnames(mat) = NULL
    # if choose to omit CI for ratio
    if(inputArg$dropRatioCi){
        mat = head(mat, 3)
    }
    ggtexttable( mat
        , theme = ttheme("blank", base_size = textsize, padding = unit(c(0.8, 0.8), "mm"))
    )
}

dCorr =
    ls(pattern = '^dat\\.', name = .GlobalEnv) %>%
    setdiff('truth') %>%
    setdiff(paste0('dat.', tolower(inputArg$ref))) %>%
    lapply(function(x){
        get(x) %>%
        mutate(Source = c(
            'dat.prev' = 'Seroprevalence'
            , 'dat.inc' = 'Seroincidence'
            , 'dat.case' = 'Case'
            , 'dat.joint' = 'Joint'
        )[x])
    }) %>%
    do.call(what = rbind) %>%
    inner_join(
        get(paste0('dat.', tolower(inputArg$ref))) %>%
            rename_at(vars(matches('[vV]al$')), function(x) paste0(x,'.ref'))
        , by = 'year'
    )
valMax = with(dCorr, max(upperVal.ref, upperVal))
valMax = foiMax = ifelse(is.null(inputArg$foiAxisMax), valMax, inputArg$foiAxisMax)

if(is.null(inputArg$foiAxisIncrement)){
    valBreaks = seq(0, round(valMax,2), by = round(valMax/4, 2))
} else {
    valBreaks = seq(0, round(valMax,2), by = inputArg$foiAxisIncrement)
}


gCorr = 
    dCorr %>%
    ggplot(aes(x = val.ref, y = val, color = Source))+
    geom_abline(slope = 1, linetype = 3, linewidth = 0.3)+
   
    geom_linerange(aes(ymin = lowerVal, ymax = upperVal), alpha = 0.15)+
    geom_linerange(aes(xmin = lowerVal.ref, xmax = upperVal.ref), alpha = 0.15)+    
    geom_point(shape = 1)+
    scale_color_manual(values = Colors$Source, guide = "none")+
    theme(
        legend.position = c(0.02,1)
        , legend.justification = c(0,1)
        , legend.background = element_rect(fill=alpha('white', 0.4))
        , axis.text.x = element_text(hjust = 0.7)
        , axis.title.y = element_text(hjust = 0.5 * valMax/foiMax)
    )+
    ylab("Estimate")+
    xlab(ifelse(inputArg$ref == 'Truth'
        , 'Truth'
        , paste0("Estimate (", inputArg$ref,")")
    ))+
    scale_x_continuous(
        breaks = valBreaks
    )+
    scale_y_continuous(
        breaks = valBreaks
    )+
    coord_cartesian(
        xlim = c(0, round(valMax + 0.01, 2))
        , ylim = c(0, round(foiMax + 0.01, 2))
    )


if(inputArg$inset){
    g =
        ggarrange(
            ggarrange(
                plotCompareValueDf(10)
                , gCorr + 
                    coord_equal(
                        ratio = 1
                        , xlim = c(0, round(valMax + 0.005, 2))
                        , ylim = c(0, round(valMax + 0.005, 2))
                    )
                , nrow = 1, ncol = 2
            )
            , gTime
            , nrow = 2
            , ncol = 1
            , heights = c(1.6,2)
        )
    ggsave(g
        , filename = outfile
        , width = 4.5
        , height = 4.5
    )  

} else {
    
    g = 
        ggarrange(gTime+
                scale_y_continuous(
                    breaks = valBreaks
                )+
                coord_cartesian(
                    ylim = c(0, round(foiMax + 0.01, 2))
                )+
                theme(axis.title.y = element_text(hjust = 0.4 * valMax/foiMax))
            , gCorr
            , nrow = 1
            , ncol = 2
            , widths = c(3,2)
            , align = "h"
        ) 
    g =
        ggdraw() +
        draw_plot(g) +
        draw_plot(plot = plotCompareValueDf()
            , x = 0.95
            , y = 1
            , hjust = 1
            , vjust = 1.2      
            , width = 0.2
            , height = 0.2
        )
    ggsave(g
        , filename = outfile
        , width = 5.5
        , height = 2
    )   
}



# CI widths
message("CI widths of FOI ....")
ls(pattern = "^dat\\.") %>%
Filter(f = function(x) !grepl('\\.plot$', x)) %>%
    lapply(function(x){
        message(x)
        get(x) %>%
        mutate(ciRange = upperVal - lowerVal) %>%
        with(summary(ciRange))
    })
