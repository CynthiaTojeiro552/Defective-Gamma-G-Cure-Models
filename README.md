# Defective Gamma-G Cure Models

This repository contains R code used in the manuscript:

**Defective Gamma-G Family for Cure Fraction Models: Novel Survival Methods with Applications to Cancer Data**

The repository provides computational routines for fitting defective Gamma-Gompertz and defective Gamma-Dagum survival models, reproducing the real-data applications, and running the Monte Carlo simulation studies reported in the manuscript.

## Repository structure

- `R/`: functions for baseline defective distributions, Gamma-G models, likelihoods, inference, simulation, and numerical routines.
- `scripts/`: scripts used to reproduce the simulation studies and the melanoma and cervical cancer applications.
- `data/`: instructions for preparing the datasets. The real datasets are not redistributed in this repository.
- `outputs/`: folder where tables, figures, and numerical results are saved.

## Simulation studies

The repository includes two simulation scripts:

- `scripts/01_simulation_gammagompertz_clean.R`: Monte Carlo simulation for the defective Gamma-Gompertz model and comparison with the defective Gompertz baseline.
- `scripts/02_simulation_gammadagum_clean.R`: Monte Carlo simulation for the defective Gamma-Dagum model and comparison with the defective Dagum baseline.

By default, the simulation scripts are set to a quick test mode:

```r
quick_run <- TRUE


## Contact

For questions about the code, please contact the corresponding author.
