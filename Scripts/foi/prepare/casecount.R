
# import age groupings
ageGroups = file.path(casedir, 'ageGroups.txt') %>% read.table
# import case count data
case = file.path(casedir, 'casecount.csv') %>% read_csv(col_type = 'd')
# import population data
pop = read_csv(popfile, col_type = 'd')

# Treat age in years as age groups.
# Expand only columns (real age) to quarterly unit
ageGroups = apply(ageGroups, 1, function(x){
    rep(x, each = 1/T.delta)
}) %>% t

# Make population annually constant
# and ages are attributed to each age increment equally
pop = lapply(colnames(pop), function(Year){
    out = lapply(pop[[Year]], matrix, ncol = 1/T.delta, nrow = 1/T.delta) %>%
        do.call(what = rbind)
    Year = as.numeric(Year)
    colnames(out) = Year + cumsum(rep(T.delta, 1/T.delta)) - T.delta
    out * T.delta
}) %>%
do.call(what = cbind)

# filter data to include only doable ages/years
# given the overlaps between case and pop data
Years = intersect(colnames(case), colnames(pop)) %>%
    as.numeric %>%
    range

# reduce population data to just those years
pop = pop[ , between(as.numeric(colnames(pop)), Years[1], Years[2])] %>%
    # and cohorts that are involved
    tail(ncol(ageGroups) + diff(Years)/T.delta - 1)

case = case[ , between(as.numeric(colnames(case)), Years[1], Years[2])]

Years = seq(Years[1] - T.pre*T.delta, Years[2] + T.delta, by = T.delta)
naVal = -1
