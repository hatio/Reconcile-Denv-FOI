
## Compare/plot estimates

All scripts in this section has `"03-plots/foi"` as the default output directory.
This can be changed via the flag `-outroot [NEW_OUTPUT PATH]`.
The first argument to all script calls are used as the output file name (within the output directory).

The sections below provides example calls to the respective scripts with input files that are available in this repository.
Please check these files for the appropriate file formats.


### Force of infection (FOI) and age-specific infection risk

Example call for temporal component of the FOI:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/compare/foi.R \
        "observed_data_jointsero/m3_Mueng.pdf" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta" \
        -fitdir.joint "02-estimates/jointsero_thresh10_20/time" \
        -ref "Case"

Example call for the age-specific infection risk (relative to the first age class):

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/compare/kappa.R \
        "observed_data_jointsero/m3_Mueng.pdf" \
        -fitdir.case "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta" \
        -fitdir.jointsero "02-estimates/jointsero_thresh10_20/time"


Both calls can take the following inputs.

|flag              |description                                        |
|------------------|---------------------------------------------------|
|`-fitdir.seroinc` | Directory of fits from seroincidence data.        |
|`-fitdir.seroprev`| Directory of fits from seroprev data.             |
|`-fitdir.case`    | Directory of fits from case count data.           |
|`-fitdir.joint`   | Directory of joint fits.                          |
|`-truthdir`       | Directory storing ground truths (for simulations).|

The reference series for the comparison can be specified using the `-ref` flag.

|flag        |description                                              |
|------------|---------------------------------------------------------|
|`-ref Inc`  | Use seroincidence-derived estimates as reference series |
|`-ref Prev` | Use seroprevalence-derived estimates as reference series|
|`-ref Case` | Use case-derived estimates as reference series          |
|`-ref Joint`| Use joint estimates as reference series                 |
|`-ref Truth`| Use simulation ground truth as reference series         |

For plots on age-specific infection risk, flags `-break.age.sero` and `-break.age.case` can be used to supply paths to
age break specification files used in the respective inferences.



### Reconstructed susceptibility

Example call to plot the comparison between susceptibility reconstructions from the FOI inferences:

    # execute code at root directory of the repository
    Rscript --vanilla Scripts/foi/compare/susceptibility.R \
        "observed_data_jointsero/m3_Mueng.pdf" \
        -S.x \
            "02-estimates/casecount/KD2021/ageMin1_ageMax99/Annual_Mueng/time/m3h.rta/S.csv" \
        -S.y \
            "02-estimates/jointsero_thresh10_20/time/S.csv" \
        -S.y.group \
            "Joint sero." \
        -lab.x "Case (extended)" \
        -lab.y "Joint serology" \
        -obs.year.start 1998 \
        -obs.year.end 2019

|flags       |descriptions                                           |
|------------|-------------------------------------------------------|
|`-S.x`      | Susceptibility in the reference inference (x-axis)    |
|`-S.y`      | Susceptibility in the comparison inferences (y-axis). Multiple can be supplied separated by a space|
|`-S.y.group`| Group names of the comparison inferences. Needs to be the same length as supplied to `S.y`|
|`-lab.x`    | x-axis label                                          |
|`-lab.y`    | y-axis label                                          |

Flags `-obs.year.start` and `-obs.year.end` can be used to specify inclusive bounds for years to be compared.
