
## Infer force of infection (FOI)

All inference scripts in this repository requires an output directory specified through the `-outbase` flag.
They also share the following default configurations for the inferences.

|flag           |default|description                                                                                                                                                                                            |
|---------------|-------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-ncpu`        | 3     | Number of cpus available (the same number of MCMC chains will be initialized).                                                                                                                        |
|`-thin_factor` | 5     | Keep every *n*-th steps/iterations in the output chains. Higher thinning leads to less posterior samples being kept which reduces correlations in the samples and the of storage needed to store them.|
|`n_iter`       | 2000  | Number of iterations to run the MCMC (warmup included).                                                                                                                                               |
|`n_warmup`     | 200   | Number of iterations treated as warmup (discarded).                                                                                                                                                   |
|`-overwriteRDS`| FALSE | Including this flag means parameters will be re-estimated even if an existing fit (stored as RDS file) is available in the output directory.                                                          |

The sections below provides example calls to the respective inference scripts with input files that are available in this repository.
Please check these files for the appropriate file formats.

**Note:** Durations of FOIs inferred are currently set to `-T.delta 1` (annual) by default. This can be changed but the behavior of the inference has not been thoroughly tested.


### FOI from longitudinal serology

Here's an example call to infer FOI from longitudinal serology.

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/seroincidence/time/fit.R \
        "01-processedData/serology/gmt" \
        -titer.threshold 10 \
        -ppos0 0 \
        -ppos1.long 0.9 \
        -ppos1.short 1 \
        -kappafile "02-estimates/jointsero_thresh10_20/time/Kappa.csv" \
        -outbase "02-estimates/seroincidence_thresh10" \
        -model.ver "set_kappaJointsero/0_0.9_1"
        
`Scripts/foi/seroincidence/time/fit.R` takes the first argument `"01-processedData/serology/gmt"` as the input directory where it expects titers of individuals given in CSV file(s) under a predefined format (see tutorial on "Simulate a cohort study" for format specification details). In brief, each file in the directory is treated as data from a distinct longitudinal study. Specific studies can be excluded from the inference using the `-studyExclude` flag. For example, `-studyExclude kps2 kps3`. By default, only bleeds collected between `-ageMin 3` and `-ageMax 30` are included. This can be changed via these respective flags.

Here, the titer threshold for seropositivity is set to 10, probability of testing positive in DENV-naives set to 0, the long-term probability of testing positive in monotypic responses set to 0.9, the probability of testing positive when infections happened during the bleeding intervals (recent) set to 1, and age-specific risk given in `"02-estimates/jointsero_thresh10_20/time/Kappa.csv"`.

The output directory here is `"02-estimates/seroincidence_thresh10/time/set_kappaJointsero/0_0.9_1"`. This is a concatenation between the output directory given to `-outbase`, the current only available inference model *time*, and the specified name of this inference attempt given via the `-model.ver` flag.



### FOI from cross-sectional serology

Example call:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/seroprevalence/${stanver}/fit.R \
        "01-processedData/serology/gmt" \
        -titer.threshold 10 \
        -ppos0 0 \
        -ppos1.long 0.9 \
        -kappafile "02-estimates/jointsero_thresh10_20/time/Kappa.csv" \
        -outbase "02-estimates/seroprevalence_thresh10" \
        -model.ver "set_kappaJointsero/0_0.9_1"

The call is highly similar to the inference from longitudinal serology except for the lack of `-ppos1.short` which is not applicable.



### FOI from case counts

Example call:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/casecount/time/fit.R \
        -model.mod "m3h.rta" \
        -popfile "00-RawData/population/KamphaengPhet.csv" \
        -casedir "01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng" \
        -outbase "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng"

The `Scripts/foi/casecount/time/fit.R` script uses regular expressions to parse the model modifications supplied to the `-model.mod` flag.
Here are examples of other model modifications that might be used.

|value          |description                   |
|---------------|------------------------------|
|`-model.mod m1h.rta`| Infection risk is uniform across ages (m1). Reporting varies in time and age (rta). *This is the default.*|
|`-model.mod m1h.rt`| Infection risk is uniform across ages (m1). Reporting varies only in time (rt).|
|`-model.mod m3h.rta`| Infection risk is heterogeneous across ages (m3). Reporting varies only in time and age (rta).|

The age breaks and year breaks for reporting (and age-specific risk) are, by default, defined by 
`"Scripts/foi/configs/ageBreaks.txt"` and `"Scripts/foi/configs/yearBreaks.txt"`.
Other files can be specified via the flags `-break.age` and `-break.year`, respectively.

For input data, the `"00-RawData/population/KamphaengPhet.csv"` provides the population size by age and year.
The `"01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng"` directory provides the case data via two files:

1. `ageGroups.txt` is a matrix that maps ages in the column names to age groups in the row names.
2. `casecount.csv` is has the years as headers and each row corresponds to counts in each of those age groups.



### FOI from multiple data types

In the manuscript, we performed **joint inference between the two serological data types** via the call:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/jointsero/time/fit.R \
        -outbase "02-estimates/jointsero_thresh10_20" \
        -indir.sero "01-processedData/serology/gmt" \
        -titer.threshold 10 20 \
        -Kappa \
        -break.age "Scripts/foi/configs/ageBreaks_cohort.txt"

All input are similar to inferences with longitudinal serology or cross-sectional serology alone except for the following flags:
- `-titer.threshold` accepting more than one threshold to help inform the titer kinetics relationship
- `-Kappa` enables estimation of age-specific infection risk (relative to the first age class)
- `-break.age` specifies the file storing age breaks between the age classes



While the manuscript did not cover **joint inference between serological and case data**, we provide here an additional inference script to perform such task.

**Important note:** As underscored in the manuscript, improper model assumptions may lead to spurious conclusions. Therefore, we advise to first perform the inferences independently prior to joint inferences.

Example call:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/joint/time/fit.R \
        -outbase "02-estimates/joint_thresh10_20" \
        -indir.sero "01-processedData/serology/gmt" \
        -casedir "01-processedData/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng" \
        -popfile "00-RawData/population/KamphaengPhet.csv" \
        -titer.threshold 10 20 \
        -Kappa \
        -break.age "Scripts/foi/configs/ageBreaks.txt"

Here, case reporting is always assumed to vary in age and time to absorb away spurious signals that may mislead the inference.
