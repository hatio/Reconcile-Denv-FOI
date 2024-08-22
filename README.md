
This repository contains the data and code used to produce results in the manuscript entitled
"Reconciling heterogeneous dengue virus infection risk estimates from different study designs."

Entry points for analyses are as follows.

1. **FOIs of observed data:** `Scripts/observed_data/__main__.sh`

2. **Simulation studies:** `Scripts/simulation/cohort/__main__.sh`

3. **Theoretical mathematical analysis:** `Scripts/theory/__main__.sh`


Directory `01-processedData` contains data that have been formatted for the analyses.
Outputs of the analyses (model estimates and summary files) are provided in the `02-estimates` directory which
are then used to produce figures and tables included the manuscript.

*Note:* Collection dates of serological data and age of the participants have been rounded to maintain anonimity of the participants (see `Scripts/observed_data/roundData.R`).


## Tutorials

The following tutorials were created to demonstrate the mechanics of the code in this repository.
We hope this would facilitate the reuse/adaption of our code for future analyses.

1. [Simulate a cohort study](Scripts/tutorials/simulate_cohort.md)

2. [Infer force of infection](Scripts/tutorials/infer.md)

3. [Compare/plot estimates](Scripts/tutorials/compare.md)
