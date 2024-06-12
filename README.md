This repository contains the data and code used to produce results in the manuscript entitled
"Reconciling heterogeneous dengue virus infection risk estimates from different study designs."

Entry points for analyses are as follows.

1. **FOIs of observed data:** `Scripts/observed_data/__main__.sh`

2. **Simulation studies:** `Scripts/simulation/cohort/__main__.sh`

3. **Theoretical mathematical analysis:** `Scripts/theory/__main__.sh`


Directory `01-processedData` contains data that have been formatted for the analyses.
Outputs of the analyses (model estimates and summary files) are provided in the `02-estimates` directory which
are then used to produce figures and tables included the manuscript.
