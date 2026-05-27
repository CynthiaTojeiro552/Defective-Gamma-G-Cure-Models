# Defective Gamma-G Cure Models

This repository contains R code used in the manuscript:

**Defective Gamma-G Family for Cure Fraction Models: Novel Survival Methods with Applications to Cancer Data**

The repository provides computational routines for fitting defective Gamma-Gompertz and defective Gamma-Dagum survival models, as well as scripts for reproducing the main real-data applications reported in the manuscript.

## Repository structure

- `R/`: functions for baseline defective distributions, Gamma-G models, likelihoods, inference, and numerical routines.
- `scripts/`: scripts used to reproduce the melanoma and cervical cancer applications.
- `data/`: instructions for preparing the datasets. The real datasets are not redistributed in this repository.
- `outputs/`: folder where tables, figures, and numerical results are saved.

## Applications

The repository includes two real-data application scripts:

- `scripts/03_application_melanoma_clean.R`: defective Gamma-Gompertz and defective Gompertz models fitted to melanoma data.
- `scripts/04_application_cervical_clean.R`: defective Gamma-Dagum and defective Dagum models fitted to cervical cancer data.

## Data availability

The real cancer registry datasets analyzed in the manuscript are not redistributed by the authors. The original data should be obtained from the official Fundação Oncocentro de São Paulo (FOSP) source. After preprocessing, users should place the required files in the `data/` folder following the instructions in `data/README_data.md`.

## Code availability

The R scripts in this repository reproduce the main likelihood-based model fitting procedures, model comparison criteria, cure-fraction estimates, estimator correlation matrices, and figures reported in the manuscript.

## Required R packages

The main R packages used are:

- `numDeriv`
- `MASS`
- `survival`
- `survminer`
- `ggplot2`
- `dplyr`
- `purrr`
- `tibble`
- `broom`
- `scales`
- `openxlsx`
- `muhaz`

## How to run

From the root folder of the repository, run:

```r
source("scripts/03_application_melanoma_clean.R")
source("scripts/04_application_cervical_clean.R")
```

The generated outputs will be saved in the `outputs/` folder.

## Contact

For questions about the code, please contact the corresponding author.
