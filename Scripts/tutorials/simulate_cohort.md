
## Simulate serological data of a cohort study

Two elements are needed to simulate serological data of a cohort study:

1. *Observation points*: Individuals being followed and time points at which their serological samples were collected

2. *Infection histories*: Ground truth of when each of the individuals were infected

Given these inputs, we can simulate titers of individuals are the *observation points* under various titer kinetics and assay noise models.

----

### Observation points

Files in the [01-processedData/serology/gmt] directory are actual data observed from each of the cohort studies in Kamphaeng Phet.
These tables can be used as observation point tables to simulate titers from equivalent cohort studies with defined infection risk, titer kinetics, and assay noise models.

We can also use the `Scripts/simulation/cohort/simulate_obs.R` script to simulate an observation points table that covers the same study periods as the actual cohort studies. Below is an example where we simulated observation points of 500,000 individuals uniformly sampled across years and ages (between 5 to 15yrs) with 2 bleeds per individual. The flag `-nogap` tells the script to not respect time gaps between the cohort studies (observations can happen in those years as well). `-ageDelta` specifies that age between bleeds of an individual are 3 months apart.

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/simulation/cohort/simulate_obs.R \
        -obs 01-processedData/serology/gmt \
        -nogap \
        -dataver 'uniform_5_15_nogap_2bleed3month_million' \
        -num_bleed 1000000 \
        -ageDelta 0.25 \
        -iFold 1

See (here)[01-processedData/simcohort/obs/uniform_5_15_nogap_2bleed3month_million/1/obs.csv] for the simulated observation points.

The observation point table contains the following headings:

- `id`: Unique identifiers of the individuals
- `dateCollect`: Dates in the standard R format (YYYY-MM-DD)
- `age`: Age of the individual at time of sample collection
- `location`: Location at which the individual resides. This is a column needed for the manuscript but is not used in the simulation studies.
- `titer`: Actual observed titers (in the case of using actual observed data as observation point table). Does not affect the simulation.
- `study`: Cohort study for which the individual was enrolled. (Study, id) pairs are used to identify unique individuals.

The flag `-iFold 1` specifies that this is realization number `1` of the observation points table.
This feature allows multiple realizations to be simulated and stored separately. 


### Infection histories

Now that we have individuals under study, the next step is to simulate infection histories of individuals using
a table of infection parameters [01-processedData/simcohort/m3/infection_param/mock] and
the observation table we generated earlier.
Note that path of the observation table(s) is supplied to the flag `-obs` as a directory which
contains the observation tables rather than the files themselves.
This is to support having observations from multiple cohorts (as is the case for Kamphaeng Phet).

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/simulation/cohort/simulate_infections.R \
        -params 01-processedData/simcohort/m3/infection_param/mock \
        -obs 01-processedData/simcohort/obs/uniform_5_15_nogap_2bleed3month_million/1 \
        -outroot 01-processedData/simcohort/m3_uniform_5_15_nogap_2bleed3month_million \
        -iFold 1

Below are the first few rows of the `01-processedData/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/infection/NA.csv` created 
from executing the code above. The `infectionAge` column are age at infections of the individual separated by semi-colons (empty if never infected).

    study,id,location,dateBirth,infectionAge
    NA,1,NA,1993-09-22,8.35995405141846;11.6935496935734
    NA,10,NA,2001-01-16,
    NA,100,NA,1991-09-16,11.2768410023729;14.769618414351
    NA,1000,NA,1997-01-12,10.6491331988941

As with simulating observation points, the flag `-iFold 1` specifies that this is realization number `1` of the infection histories table.
This feature allows multiple realizations to be simulated and stored separately.


### Simulate titers

The code below simulates observed titers at observation points supplied via the directory of the observation tables and
infection histories under the subdirectory `infections` in the directory specified via the `-outroot` flag.
The `-scenario` flag provides a name for the titer kinetics and assay noise configurations used to simulate the titers.

Available configurations:

- `-Omega`: Permanent rise in titers upon 1st DENV infection
- `-Gamma`: Temporary rise in titers upon 1st DENV infection
- `-Delta`: Exponential decay rate of the temporary rise (per year)
- `-Omega0_rel`: Cross-reactive titers in DENV-naive individuals as multiples of `Omega`
- `-Sigma`: Standard deviation of the normally-distributed assay noise

Note that these titers are on adjusted log2 scales: `log2(linear titer/10) + 1`.

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/simulation/cohort/simulate_titers.R \
        -obs 01-processedData/simcohort/obs/uniform_5_15_nogap_2bleed3month_million/1 \
        -outroot 01-processedData/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1 \
        -scenario noisepos_steady \
        -Omega0_rel 0.2 -Delta 0

For configurations not specified, (default values)[Scripts/simulation/cohort/simulate_titers.R] applies.

The above example generates an observed titer table stored in
`01-processedData/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noisepos_steady`.

Here are the first few rows of what this file might look like:

    dateCollect,id,age,titer,study,Titer.true
    1999-02-09,1,5.379799754184205,6.383099278491924,NA,6.012346244490888
    1999-05-11,1,5.629799754184205,5.35861099905044,NA,6.012346244490888
    2007-02-15,10,6.079542569292244,4.4744763497580315,NA,6.012346244490888
    
The columns `dateCollect`, `id`, `age`, and `study` are taken from the input observation table(s).
`titer` is the simulated observed titer while `Titer.true` is the true underlying titer of the individual at that point
based on the assumed titer kinetics, before applying assay noise.

