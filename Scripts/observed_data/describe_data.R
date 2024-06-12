
library(tidyverse); theme_set(theme_classic())
library(ggpubr)


source("Scripts/configs/Aesthetics.R")

    #   Import data
    #   ...........

casedir = c(
    'All'     = "01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_All"
    , 'Mueng' = "01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng"
)
popfile = "01-processedData/population/KamphaengPhet.csv"
indir = "01-processedData/serology/gmt"

# ... case data

    # import case count data
    case = 
        casedir %>%
        lapply(function(casedir){
            file.path(casedir, 'casecount.csv') %>% read_csv(col_type = 'd')
        }) %>%
        lapply(function(case){
            case %>%
            mutate(Age = 1:nrow(case) - 1) %>%
            pivot_longer(-Age, names_to = "Year", values_to = "count") %>%
            filter(count != -1) %>%
            mutate(Year = as.integer(Year))          
        }) 
    # import population data
    pop = read_csv(popfile)

    pop = pop %>%
        mutate(
            cohort = as.integer(tail(colnames(pop),1)) - nrow(pop):1 + 1
        ) %>%
        pivot_longer(-cohort, names_to = "Year", values_to = "pop") %>%
        mutate(Year = as.integer(Year), Age = Year - cohort) %>%
        select(-cohort)

    case = 
        case$All %>%
        inner_join(case$Mueng
            , by = c('Year', 'Age')
            , suffix = c('.all', '.mueng')
        ) %>%
        inner_join(pop
            , by = c('Year', 'Age')
        )

# ... serology

    dat = 
        list.files(indir, full.names = T) %>%
        lapply(read_csv) %>%
        do.call(what = rbind)

    dat = dat %>%
        filter(!is.na(titer)) %>%
        mutate(
            # convert dateCollect to fractional years
            year = as.numeric(format(dateCollect,'%Y')) +
                as.numeric(dateCollect - as.Date(paste0((format(dateCollect,'%Y')),'-01-01')))/365.25
        ) %>%
        select(study, id, age, year, titer)



    #   Describe data
    #   .............

# %sero+ at age 9
dat %>%
    filter(between(age, 9, 10)) %>%
    group_by(study,id) %>%
    filter(age == max(age)) %>%
    group_by(study) %>%
    summarize(
        num = n()
        , ppos = round(sum(titer >= 10)*100 / n())
    )

# %sero+ after age 30 in KFCS
dat %>%
    filter(age >= 31, study == 'kfcs') %>%
    group_by(study,id) %>%
    filter(age == min(age)) %>%
    group_by(study) %>%
    summarize(
        ppos.10 = round(sum(titer >= 10)*100 / n())
        , ppos.20 = round(sum(titer >= 20)*100 / n())
    )

# number of cases (mueng vs all)
case %>%
    summarize(
        sum(count.all)
        , sum(count.mueng)
        , prop.mueng = round(sum(count.mueng)/sum(count.all) * 100)
    )

# case per 1000
case %>%
    group_by(Year) %>%
    summarize(
        casePer1000 = round(sum(count.all) / sum(pop) * 1000,1)
    ) %>%
    arrange(casePer1000) %>%
    as.data.frame

# mean age of cases in Mueng
case %>%
    group_by(Year) %>%
    summarize(
        meanAge = round(sum(Age * count.mueng) / sum(count.mueng),1)
    ) %>%
    as.data.frame


    #   Plot data
    #   .........

titer.threshold = 10
# yearLimits = c(1993.5,2019.5)
yearLimits = c(1994,2020)
ageMin = 3; ageMax = 30

gTs.sero =
    dat %>%
    mutate(Year = floor(year)) %>%
    group_by(Year) %>%
    summarize(count = n()) %>%
    mutate(Year = Year - 0.5) %>%
    ggplot(aes(x = Year, y = count))+
    geom_col(aes(fill = ""), width = 1, position = "stack")+
    scale_fill_manual("", values = Colors$Source[['Seroprevalence']])+
    scale_x_continuous(expand = c(0,0), limits = yearLimits)+
    scale_y_continuous("Bleeds", expand = c(0,0))

gSero = 
    dat %>%
    filter(age < (ageMax+1)) %>%
    mutate(pos = titer >= titer.threshold) %>%
    mutate(Age = floor(age), Year = floor(year) + 0.5) %>%
    group_by(Age, Year) %>%
    # remove ages with few data points (unreliable proportions)
    filter(n() >= 10) %>%
    summarize(Prop = mean(pos)) %>%
    ggplot(aes(x = Year, y = Age))+
    geom_tile(aes(fill = Prop))+
    scale_fill_distiller("% Seropositive"
        , palette = "YlOrBr", direction = 1
        , labels = scales::percent_format(scale = 100)
        , breaks = c(0,0.5,1)
        , guide = guide_colourbar(
            direction = "horizontal"
            , title.position = "top"
            , barheight = rel(0.5)
        )
    )+
    scale_x_continuous(expand = c(0,0), limits = yearLimits)+
    scale_y_continuous(expand = c(0,0))

gCase = 
    case %>%
    mutate(Year = Year - 0.5) %>%
    ggplot(aes(x = Year, Age))+
    geom_tile(aes(fill = count.all/pop * 10^3))+
    scale_fill_continuous("Cases per 1000"
        , low = "#eae9f2", high = Colors$Source[['Case']]
        , guide = guide_colourbar(
            direction = "horizontal"
            , title.position = "top"
            , barheight = rel(0.5)
        )
    )+
    scale_x_continuous(expand = c(0,0), limits = yearLimits)+
    scale_y_continuous(expand = c(0,0))
    
gTs.case = 
    case %>%
    group_by(Year) %>%
    summarize(count = sum(count.all), pop = sum(pop)) %>%
    mutate(Year = Year - 0.5) %>%
    ggplot(aes(x = Year, y = count/pop * 10^3))+
    geom_col(aes(fill = ""), width = 1)+
    scale_fill_manual("", values = Colors$Source[['Case']])+
    scale_x_continuous(expand = c(0,0), limits = yearLimits)+
    scale_y_continuous("Case\nper 1000", expand = c(0,0))



    #   Table S1
    #   ........

# age, year ranges
dat %>%
    filter(between(age, ageMin, ageMax)) %>%
    select(study, year, age) %>%
    group_by(study) %>%
    summarize_at(vars(everything()), range) %>%
    as.data.frame

# number of bleeds, individuals
dat %>%
    filter(between(age, ageMin, ageMax)) %>%
    group_by(study, id) %>%
    summarize(numBleed = n()) %>%
    group_by(study) %>%
    summarize(
        numBleed = sum(numBleed)
        , numIndiv = n()
    )
# number of family clusters (KFCS)
dat %>%
    filter(between(age, ageMin, ageMax), study == 'kfcs') %>%
    mutate(clusterId = str_extract(id, '^[0-9]{3}KF')) %>%
    with(unique(clusterId)) %>%
    length
    



    
    #   Plot map
    #   ........

spacedir = "00-RawData/spatial"

# # download shapefile
# dir.create(spacedir, recursive = T)
# geodata::gadm(
#     country = 'THA', level = 3
#     , path = spacedir
# )

# import Thai map
map = 
    file.path(spacedir, "gadm36_THA_3_sp.rds") %>% 
    readRDS %>%
    sf::st_as_sf() %>%
    # reduce to just KPP
    filter(NAME_1 == "Kamphaeng Phet") %>%
    # merge subdistricts to form districts
    group_by(NAME_2) %>%
    summarize(geometry = sf::st_union(geometry))

# coordinate of KPPH
kpph.coord = data.frame(lat = 16.462651, long = 99.525646)

# study area
spaceAlpha = c(
    "Khanu Woralaksaburi" = 0.25
    , "Mueang Kamphaeng Phet" = 1    
)

gMap =
    map %>%
    ggplot()+
    # fill study districts
    geom_sf(data = map %>% filter(NAME_2 %in% names(spaceAlpha))
        , aes(alpha = NAME_2)
        , fill = Colors$Source[['Seroprevalence']]
    )+
    geom_sf(color = "black", fill = NA, linewidth = 0.3)+
    geom_point(data = kpph.coord, aes(x = long, y = lat), size = 7, color = Colors$Source[['Case']])+
    scale_alpha_manual(values = spaceAlpha)+
    # use WGS84 projection
    coord_sf(crs = sf::st_crs(4326))+
    theme_void()+
    theme(legend.position = "none")+
    ggsn::scalebar(data = map, transform = T, model = "WGS84"
        , dist = 20, dist_unit = "km"
        , location = 'topright'
        , st.size = 4
        , border.size = 0.5
    )




    #   Combine plots into panel
    #   .............

    outdir = "03-plots/manuscript"
    dir.create(outdir, recursive = T)


# Main text figure
g =
    ggarrange(
        gMap
        , ggarrange(
            gTs.sero %+% theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")
            , gSero %+% theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = c(0.02, 1.05), legend.justification = c(0,1), legend.background = element_rect(fill = NA))
            , gTs.case %+% theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")
            , gCase %+% theme(axis.text.x = element_text(hjust = 0.8), legend.position = c(0.01, 1.05), legend.justification = c(0,1), legend.background = element_rect(fill = NA))

            , nrow = 4
            , ncol = 1
            , heights = c(0.6, 1.5,0.6,1.5)
            , align = "v"
            , vjust = 0.5
            , hjust = 0
            , labels = c("b", "", "c", "")
        )
        , nrow = 2
        , ncol = 1
        , heights = c(2,4)
        , labels = c('a', '')
        , hjust = 0
    )
ggsave(g
    , filename = file.path(outdir, "Figure_1.pdf")
    , width = 3
    , height = 7
)

