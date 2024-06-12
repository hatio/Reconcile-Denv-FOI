
fitSeroprevalence_observed(){
    cmd="Rscript --vanilla Scripts/foi/seroprevalence/${stanver}/fit.R \
        "${datadir}" \
        -outbase "02-estimates/seroprevalence_thresh${titer_thresh}" \
        -titer.threshold "${titer_thresh}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        "
    echo '---------------'
    sbatch \
        --job-name="seroprev_${stanver}:observed" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE}
    echo Command: $cmd
    echo '---------------'
}

fitSeroincidence_observed(){
    cmd="Rscript --vanilla Scripts/foi/seroincidence/${stanver}/fit.R \
        "${datadir}" \
        -outbase "02-estimates/seroincidence_thresh${titer_thresh}" \
        -titer.threshold "${titer_thresh}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        "
    echo '---------------'
    sbatch \
        --job-name="seroinc_${stanver}:observed" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE}
    echo Command: $cmd
    echo '---------------'
}

fitCasecount_observed(){
    cmd="Rscript --vanilla Scripts/foi/casecount/${stanver}/fit.R \
        -casedir "${casedatadir}" \
        -popfile "00-RawData/population/KamphaengPhet.csv" \
        -outbase "02-estimates${casedatadir##01-processedData}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        "
    echo '---------------'
    sbatch \
        --job-name="case_${stanver}:observed" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE}
    echo Command: $cmd
    echo '---------------'
}

fitJointSeroModel_observed(){
    cmd="Rscript --vanilla Scripts/foi/jointsero/${stanver}/fit.R \
        -outbase "02-estimates/jointsero_thresh${titer_thresh/ /_}" \
        -indir.sero "${datadir}" \
        -casedir "${casedatadir}" \
        -popfile "00-RawData/population/KamphaengPhet.csv" \
        -titer.threshold "${titer_thresh}" \
        -ncpu ${NUM_CHAINS} \
        ${FLAG_RSCRIPT} \
        "
    jobId=$( sbatch \
        --job-name="joint_${stanver}:observed" \
        --export=modules_to_load=${MODULES_TO_LOAD},cmd="${cmd}" \
        --mem-per-cpu=3gb \
        --cpus-per-task=${NUM_CHAINS} \
        ${FLAG_SLURM} \
        ${@} \
        ${SBATCH_FILE} | grep -oE '[0-9]+$' )
    echo '---------------'
    echo Command: $cmd
    echo Submitted job: $jobId
    echo '---------------'
}


# function to combine estimates from all folds
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


# Path to datasets
casedatadir="01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng"
datadir="01-processedData/serology/gmt"






    #   Describe data
    #   .............
    
Rscript --vanilla Scripts/observed_data/describe_data.R





    #   Fit standard models
    #   ...................

# .... Serology
for titer_thresh in 10 20 ; do
    # fit FOI from seroprevalence
    stanver="time"
    NUM_CHAINS=5
    slurmflag_seroprev="--array=1-3"
    FLAG_RSCRIPT="-model.ver simple"
    fitSeroprevalence_observed ${slurmflag_seroprev}

    # fit FOI from seroincidence
    stanver="time"
    NUM_CHAINS=5
    FLAG_RSCRIPT="-model.ver simple"
    fitSeroincidence_observed
done


# ..... Case data
stanver="time"
NUM_CHAINS=5
FLAG_RSCRIPT="-model.mod m1h.rta"
fitCasecount_observed





    #   Fit extended models
    #   ...................

# .... Joint serology (seroprevalence, seroincidence)
FLAG_MODEL_MOD="-Kappa"
stanver="time"
titer_thresh="10 20"
NUM_CHAINS=5
flags="--array=1-3"; FLAG_RSCRIPT="${FLAG_MODEL_MOD}"
fitJointSeroModel_observed ${flags}


# .... Case data
stanver="time"
NUM_CHAINS=5
FLAG_RSCRIPT="-model.mod m3h.rta"
fitCasecount_observed


    #   Demonstrate sensitivity of serology-derived FOI
    #   to assumed test positivity probabilities
    #   : Use Kappa(a) inferred from case data
    #   ......................................


# Kappa(a) from case data
Kappa_ver="Case"
KappaFile="02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta/Kappa.csv"

# Kappa(a) from joint serology
Kappa_ver="Jointsero"
KappaFile="02-estimates/jointsero_thresh10_20/time/Kappa.csv"


# .... fit FOI from seroprevalence
stanver="time"
NUM_CHAINS=5
slurmflag_seroprev="--array=1-3"

    # Scan through combinations of test positive probabilities
    ppos_ver="set"
    for titer_thresh in 10 20 ; do
    for ppos0 in $( seq -w 0 0.02 0.1); do
    for ppos1_long in $( seq -w 0.9 0.02 1); do

        FLAG_RSCRIPT="-model.ver ${ppos_ver}_kappa${Kappa_ver} -ppos0 $ppos0 -ppos1.long $ppos1_long -kappafile ${KappaFile}"
        fitSeroprevalence_observed ${slurmflag_seroprev}

    done
    done
    done


# .... fit FOI from seroincidence
stanver="time"
NUM_CHAINS=5

    # Scan through combinations of test positive probabilities
    ppos_ver="set"
    for titer_thresh in 10 20 ; do
    for ppos0 in $( seq -w 0 0.02 0.1); do
    for ppos1_long in $( seq -w 0.9 0.02 1); do

        FLAG_RSCRIPT="-model.ver ${ppos_ver}_kappa${Kappa_ver} -ppos0 $ppos0 -ppos1.long $ppos1_long -kappafile ${KappaFile}"
        fitSeroincidence_observed

    done
    done
    done



    #   Combine estimates across folds
    #   ..............................

# seroprevalence
for fitdir in $( find 02-estimates -type f | \
    grep 'fit.RDS$' | \
    grep 'seroprevalence_thresh' | \
    sed -E 's/\/[0-9]+\/fit\.RDS$//g' | \
    uniq ) ; do
    
    combineFold $fitdir
done

# jointsero
for fitdir in $( find 02-estimates -type f | \
    grep 'fit.RDS$' | \
    grep 'jointsero_thresh' | \
    sed -E 's/\/[0-9]+\/fit\.RDS$//g' | \
    uniq ) ; do
    
    combineFold $fitdir
done





    #   Plot results from standard model
    #   ................................

case_area="Mueng"
titer_thresh=10
titer_thresh=20

    # Annual FOI
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_thresh${titer_thresh}/m1_${case_area}.pdf" \
        -fitdir.seroinc "02-estimates/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta" \
        -ref "Case" \
        -salje
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_thresh${titer_thresh}/m1_${case_area}_dropCi.pdf" \
        -fitdir.seroinc "02-estimates/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta" \
        -ref "Case" \
        -dropRatioCi \
        -salje
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_thresh${titer_thresh}/m1_${case_area}_refPrev.pdf" \
        -fitdir.seroinc "02-estimates/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta" \
        -ref "Prev"

    # Reconstructed susceptibility (vs case)
    Rscript --vanilla Scripts/foi/compare/susceptibility.R \
        "observed_data_thresh${titer_thresh}/m1_${case_area}_vs_case.pdf" \
        -S.x \
            "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta/S.csv" \
        -S.y \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/1/S.csv" \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/2/S.csv" \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/3/S.csv" \
            "02-estimates/seroincidence_thresh${titer_thresh}/time/simple/S.csv" \
        -S.y.group \
            "Seroprevalence" \
            "Seroprevalence" \
            "Seroprevalence" \
            "Seroincidence" \
        -lab.x "Case-derived" \
        -lab.y "Serology-derived"


    # Annual FOI
    # (only lab-confirmed cases for case-derived FOIs)
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_thresh${titer_thresh}/labconfirmed_m1_${case_area}.pdf" \
        -fitdir.seroinc "02-estimates/seroincidence_thresh${titer_thresh}/time/simple" \
        -fitdir.seroprev "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple" \
        -fitdir.case "02-estimates/casecount/KD2021_dengue/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta" \
        -ref "Case"

    # Reconstructed susceptibility (vs case)
    # (only lab-confirmed cases for case-derived FOIs)
    Rscript --vanilla Scripts/foi/compare/susceptibility.R \
        "observed_data_thresh${titer_thresh}/lab_confirmed_m1_${case_area}_vs_case.pdf" \
        -S.x \
            "02-estimates/casecount/KD2021_dengue/ageMin1_ageMax99/Annual_${case_area}/time/m1h.rta/S.csv" \
        -S.y \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/1/S.csv" \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/2/S.csv" \
            "02-estimates/seroprevalence_thresh${titer_thresh}/time/simple/3/S.csv" \
            "02-estimates/seroincidence_thresh${titer_thresh}/time/simple/S.csv" \
        -S.y.group \
            "Seroprevalence" \
            "Seroprevalence" \
            "Seroprevalence" \
            "Seroincidence" \
        -lab.x "Case-derived" \
        -lab.y "Serology-derived"


    #   Plot results from extended models
    #   .................................

case_area="Mueng"

    # FOI
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_jointsero/m3_${case_area}.pdf" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m3h.rta" \
        -fitdir.joint "02-estimates/jointsero_thresh10_20/time" \
        -ref "Case"

    # Kappa(a)
    Rscript --vanilla Scripts/foi/compare/kappa.R \
        "observed_data_jointsero/m3_${case_area}.pdf" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m3h.rta" \
        -fitdir.jointsero "02-estimates/jointsero_thresh10_20/time"

    # Reconstructed susceptibility
    Rscript --vanilla Scripts/foi/compare/susceptibility.R \
        "observed_data_jointsero/m3_${case_area}.pdf" \
        -S.x \
            "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_${case_area}/time/m3h.rta/S.csv" \
        -S.y \
            "02-estimates/jointsero_thresh10_20/time/1/S.csv" \
            "02-estimates/jointsero_thresh10_20/time/2/S.csv" \
            "02-estimates/jointsero_thresh10_20/time/3/S.csv" \
        -S.y.group \
            "Joint sero." \
            "Joint sero." \
            "Joint sero." \
        -lab.x "Case (extended)" \
        -lab.y "Joint serology"



    #   Plot sensitivity of serology-derived FOI
    #   ........................................
  
# titer thresh = 10  
Rscript --vanilla Scripts/foi/compare/ppos_foi.R \
    "observed_data_thresh10/Mueng_kappaJointsero.pdf" \
    -fitdir.seroinc "02-estimates/seroincidence_thresh10/time/set_kappaJointsero" \
    -fitdir.seroprev "02-estimates/seroprevalence_thresh10/time/set_kappaJointsero" \
    -fitdir.ref "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta" \
    -ppos0 0.07203656416889753 \
    -ppos1.long 0.9996810774944899

# titer thresh = 20
Rscript --vanilla Scripts/foi/compare/ppos_foi.R \
    "observed_data_thresh20/Mueng_kappaJointsero.pdf" \
    -fitdir.seroinc "02-estimates/seroincidence_thresh20/time/set_kappaJointsero" \
    -fitdir.seroprev "02-estimates/seroprevalence_thresh20/time/set_kappaJointsero" \
    -fitdir.ref "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta" \
    -ppos0 3.061701422826246e-4 \
    -ppos1.long 0.9256323037252381


    #   Format estimates for supplementary tables
    #   .........................................
    
Rscript --vanilla Scripts/foi/format/foi.R "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta"
Rscript --vanilla Scripts/foi/format/foi.R "02-estimates/jointsero_thresh10_20/time"
