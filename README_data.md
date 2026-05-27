# Data instructions

The real cancer registry datasets used in the manuscript are not redistributed in this repository.

Users interested in reproducing the real-data applications should obtain the original data from the official Fundação Oncocentro de São Paulo (FOSP) source and prepare the datasets with the variable names expected by the scripts.

## Melanoma application

The script `scripts/03_application_melanoma_clean.R` expects a file named:

`data/melanoma.txt`

with at least the following columns:

- `tempo_anos`: observed survival or censoring time in years;
- `status`: event indicator, with 1 for event and 0 for censoring;
- `EC_cat`: clinical stage group indicator.

## Cervical cancer application

The script `scripts/04_application_cervical_clean.R` expects a file named:

`data/cervical.txt`

with at least the following columns:

- `tempo`: observed survival or censoring time;
- `status`: event indicator, with 1 for event and 0 for censoring;
- `ECGRUP`: clinical stage group.

In the cervical cancer application, stages I and II are grouped as `I_II`, while stages III and IV are grouped as `III_IV` within the analysis script.

## Important note

The files `melanoma.txt` and `cervical.txt` are not included in the public repository and should not be uploaded if data redistribution is not permitted.
