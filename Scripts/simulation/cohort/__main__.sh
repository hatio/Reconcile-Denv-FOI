
    #   Prepare observation points (age, year of bleeds)
    #   ..........................

# uniform sampling across years and age (between 5 to 15), without gap years between cohorts, million bleeds, 2 bleeds per individual
Rscript --vanilla Scripts/simulation/cohort/simulate_obs.R -nogap -dataver 'uniform_5_15_nogap_2bleed3month_million' -num_bleed 1000000 -ageDelta 0.25

# uniform sampling across years and age (between 5 to 15), without gap years between cohorts, million bleeds, 4 bleeds per individual
Rscript --vanilla Scripts/simulation/cohort/simulate_obs.R -nogap -dataver 'uniform_5_15_nogap_4bleed3month_million' -num_bleed 1000000 -num_bleed_per_indiv 4 -ageDelta 0.25

# uniform sampling across years and age (between 5 to 15), without gap years between cohorts, 500k bleeds, 4 bleeds per individual
Rscript --vanilla Scripts/simulation/cohort/simulate_obs.R -nogap -dataver 'uniform_5_15_nogap_4bleed3month_500k' -num_bleed 500000 -num_bleed_per_indiv 4 -ageDelta 0.25

# copy  actual observation points
mkdir -p "01-processedData/simcohort/obs/actual/1"
for f in $(ls 01-processedData/serology/gmt ); do
    cp -s ${PWD}/01-processedData/serology/gmt/$f "${PWD}/01-processedData/simcohort/obs/actual/1/$f"
done




    #   (Function to) Simulate infection histories
    #   ..........................................

simulateInfections(){
    if [ "$foiModel" == "m1" ]; then
        foiDir="01-processedData/simcohort/m1/infection_param/mock"
    elif [ "$foiModel" == "m3" ]; then
        foiDir="01-processedData/simcohort/m3/infection_param/mock"
    elif [ "$foiModel" == "m0" ]; then
        foiDir="01-processedData/simcohort/m0/infection_param/mock"
    else
        echo "Invalid FOI model: ${foiModel}"
        return 0
    fi
    cmd="Rscript --vanilla Scripts/simulation/cohort/simulate_infections.R \
        -params ${foiDir} \
        -obs 01-processedData/simcohort/obs/${obsver}/1 \
        -outroot 01-processedData/simcohort/${foiModel}_${obsver} \
        -iFold $iFold
        "
    jobId=$( sbatch \
        --job-name="infection:${foiModel}_${obsver}" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --cpus-per-task=1 \
        --mem-per-cpu=3gb \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE} | grep -oE '[0-9]+$' )
    echo '---------------'
    echo Submitted job: $jobId
    echo Command: $cmd
    echo '---------------'
    echo ''
}


    #   (Function to) Simulate titers
    #   .............................

simulateTiters(){
    if [ "$scenario" == "perfect_steady" ]; then
        local flag_scenario="-Delta 0 -Sigma 0 -Omega0_rel 0"
    elif [ "$scenario" == "perfect_wane" ]; then
        local flag_scenario="-Sigma 0 -Omega0_rel 0"
    elif [ "$scenario" == "noise_steady" ]; then
        local flag_scenario="-Delta 0"
    elif [ "$scenario" == "noise_wane" ]; then
        local flag_scenario=""
    elif [ "$scenario" == "noisepos_steady" ]; then
        local flag_scenario="-Delta 0 -Omega0_rel 0"
    elif [ "$scenario" == "noisepos_wane" ]; then
        local flag_scenario="-Omega0_rel 0"
    else
        echo "Invalid titer scenario: ${scenario}"
        return 0
    fi
    cmd="Rscript --vanilla Scripts/simulation/cohort/simulate_titers.R \
        -obs 01-processedData/simcohort/obs/${obsver}/1 \
        -outroot 01-processedData/simcohort/${foiModel}_${obsver}/${iFold} \
        -scenario ${scenario} \
        ${flag_scenario} \
        "
    jobId=$( sbatch \
        --job-name="titer:${foiModel}_${obsver}:${scenario}" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --cpus-per-task=1 \
        --mem-per-cpu=3gb \
        ${FLAG_SLURM} \
        --dependency=afterany:${jobId} \
        ${@} \
        ${SBATCH_FILE} | grep -oE '[0-9]+$' )
    echo '---------------'
    echo Submitted job: $jobId
    echo Command: $cmd
    echo '---------------'
    echo ''
}


    #   (Functions to) infer FOI
    #   ........................

fitSeroprevalence_sim(){
    cmd="Rscript --vanilla Scripts/foi/seroprevalence/${stanver}/fit.R \
        "01-processedData/simcohort/${foiModel}_${obsver}/${iFold}/${scenario}" \
        -outbase "02-estimates/simcohort/${foiModel}_${obsver}/${iFold}/${scenario}/seroprevalence_thresh${titer_thresh}" \
        -titer.threshold "${titer_thresh}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        ${FLAG_OVERWRITE} \
        "
    echo '---------------'
    sbatch \
        --job-name="seroprev_${stanver}:${foiModel}_${obsver}:${scenario}" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE}
    echo Command: $cmd
    echo '---------------'
}

fitSeroincidence_sim(){
    cmd="Rscript --vanilla Scripts/foi/seroincidence/${stanver}/fit.R \
        "01-processedData/simcohort/${foiModel}_${obsver}/${iFold}/${scenario}" \
        -outbase "02-estimates/simcohort/${foiModel}_${obsver}/${iFold}/${scenario}/seroincidence_thresh${titer_thresh}" \
        -titer.threshold "${titer_thresh}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        ${FLAG_OVERWRITE} \
        "
    echo '---------------'
    sbatch \
        --job-name="seroinc_${stanver}:${foiModel}_${obsver}:${scenario}" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE}
    echo Command: $cmd
    echo '---------------'
}



    #   (Functions to) combine estimates from all folds
    #   ........................

combineFold(){
    cmd="Rscript --vanilla Scripts/foi/summarize/combineFold.R ${1}"
    sbatch \
        --job-name="combine" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=10gb \
        ${SBATCH_FILE}
    echo Combine fold command: $cmd
}



# FLAGS related to SLURM
MODULES_TO_LOAD="R/4.1"
FLAG_SLURM="--time=96:00:00 --mem-per-cpu=15gb"
SBATCH_FILE="Scripts/configs/universal.sbatch"



    #   Simulate data
    #   .............

# simulate infections
tail -qn +2 01-processedData/simcohort/thresh*.csv | \
cut -d"," -f2-4 | sort -r | uniq | \
while IFS=, read -r iFold obsver foiModel
do
    simulateInfections
done


# simulate titers
tail -qn +2 01-processedData/simcohort/thresh*.csv | \
cut -d"," -f2-5 | sort -r | uniq | \
while IFS=, read -r iFold obsver foiModel scenario
do
    jobId=1
    simulateTiters
done








    #   Fit models
    #   ..........

FLAG_OVERWRITE=""
NUM_CHAINS=5

# fit models: estimate foi
stanver="time"
control_file="thresh10"
control_file="thresh20"
control_file="thresh10_actual"
control_file="thresh20_actual"



    # seroincidence
    sed '1d' 01-processedData/simcohort/${control_file}.csv | \
    while IFS=, read -r titer_thresh iFold obsver foiModel scenario rflag_seroinc rflag_seroprev slurmflag_seroprev analysis_set
    do
        FLAG_RSCRIPT=${rflag_seroinc}
        if [[ "$FLAG_RSCRIPT" == "" ]]; then
            continue
        fi
        fitSeroincidence_sim
    done

    # seroprevalence
    sed '1d' 01-processedData/simcohort/${control_file}.csv | \
    while IFS=, read -r titer_thresh iFold obsver foiModel scenario rflag_seroinc rflag_seroprev slurmflag_seroprev analysis_set
    do
        FLAG_RSCRIPT=${rflag_seroprev}
        if [[ "$FLAG_RSCRIPT" == "" ]]; then
            continue
        fi
        fitSeroprevalence_sim ${slurmflag_seroprev}
    done
    


    #   Combine estimates across folds
    #   ..............................

# seroprevalence
for fitdir in $( find 02-estimates_new -type f | \
    grep 'fit.RDS$' | \
    grep 'seroprevalence_thresh' | \
    sed -E 's/\/[0-9]+\/fit\.RDS$//g' ) ; do
    
    combineFold $fitdir
done






    #   Plot FOI estimates (new narrative)
    #   ..................

titer_thresh=10
titer_thresh=20

# Ideal scenario
Rscript --vanilla Scripts/foi/compare/foi.R \
    "simcohort_thresh${titer_thresh}/high_power/ideal.pdf" \
    -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/simple" \
    -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/simple"

# Effects of violating ONE assumption

    # Waning monotypic titers
    # So long as exposure status of individuals could be perfectly ascertained, it still works
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_1/perfect_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_wane/seroprevalence_thresh${titer_thresh}/time/simple"

    # When assay is imperfect (noisy)
    # - seroprev is still okay
    # - seroinc becomes inflated
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_1/noisepos_steady.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_steady/seroprevalence_thresh${titer_thresh}/time/simple"

    # Effects of neglecting non-uniform risk in age
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_1/kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/simple" \
        -truthdir "01-processedData/simcohort/m3/infection_param/mock"

# Effects of violating MULTIPLE assumptions
# : seroinc inflation exacerbated 

    # when assay is noisy & waning monotypic titers
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_multi/noisepos_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_wane/seroprevalence_thresh${titer_thresh}/time/simple"

    # when assay is noisy & cxr titers
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_multi/noise_steady.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_steady/seroprevalence_thresh${titer_thresh}/time/simple"

    # when assay is noisy & waning monotypic titers & cxr
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_multi/noise_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/simple"

    # when assay is noisy & waning monotypic titers & cxr & non-uniform risk in age
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/violate_multi/noise_wane_kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/simple"

# Re-plot these effects for main text figure
titer_thresh=10

    # Ideal scenario
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/panel_c_high_power_ideal.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/simple" \
        -foiAxisMax 0.42

    # when assay is noisy & waning monotypic titers & cxr
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/panel_d_high_power_noise_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/simple" \
        -foiAxisMax 0.42
    
    # Effects of neglecting non-uniform risk in age
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/panel_e_high_power_kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/simple" \
        -truthdir "01-processedData/simcohort/m3/infection_param/mock" \
        -foiAxisMax 0.42

    # when assay is noisy & waning monotypic titers & cxr & non-uniform risk in age
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/panel_f_high_power_noise_wane_kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/simple" \
        -foiAxisMax 0.42


# Efficiency of correcting violations

    # assay noise alone
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/relax/noisepos_steady.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_steady/seroincidence_thresh${titer_thresh}/time/noisepos" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_steady/seroprevalence_thresh${titer_thresh}/time/noisepos"

    # assay noise & waning monotypic titers
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/relax/noisepos_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_wane/seroincidence_thresh${titer_thresh}/time/noisepos_wane" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noisepos_wane/seroprevalence_thresh${titer_thresh}/time/noisepos_wane"

    # assay noise & waning monotypic titers & cxr titers
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/relax/noise_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/noise_wane" \
        -fitdir.seroprev "02-estimates/simcohort/m1_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/noise_wane"

    # assay noise & waning monotypic titers & cxr titers & non-uniform risk in age
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/relax/noise_wane_kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroincidence_thresh${titer_thresh}/time/noise_wane_kappa" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/noise_wane_kappa"
    # the correction inefficiency when risk is non-uniform in age is similar to when risk is uniform in age
    
    # non-uniform risk in age alone
    # can be efficiently corrected for
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/high_power/relax/kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/kappa" \
        -fitdir.seroprev "02-estimates/simcohort/m3_uniform_5_15_nogap_2bleed3month_million/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/kappa" \
        -truthdir "01-processedData/simcohort/m3/infection_param/mock"


# Reducing power to match cohort studies
titer_thresh=10
titer_thresh=20

    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/actual_power/relax/relax_noise_wane.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_actual/1/noise_wane/seroincidence_thresh${titer_thresh}/time/noise_wane" \
        -fitdir.seroprev "02-estimates/simcohort/m1_actual/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/noise_wane" \
        -truthdir "01-processedData/simcohort/m1/infection_param/mock"

    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/actual_power/relax/relax_noise_wane_kappa.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m3_actual/1/noise_wane/seroincidence_thresh${titer_thresh}/time/noise_wane_kappa" \
        -fitdir.seroprev "02-estimates/simcohort/m3_actual/1/noise_wane/seroprevalence_thresh${titer_thresh}/time/noise_wane_kappa" \
        -truthdir "01-processedData/simcohort/m3/infection_param/mock"

    # In fact, uncertainty is still very high (especially in seroprevalence) even when no assumption violated
    # (perfect assay, no waning titers, uniform risk in age)
    # : uncertainty purely due to reduced power
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "simcohort_thresh${titer_thresh}/actual_power/ideal.pdf" \
        -fitdir.seroinc "02-estimates/simcohort/m1_actual/1/perfect_steady/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/simcohort/m1_actual/1/perfect_steady/seroprevalence_thresh${titer_thresh}/time/simple"
