# anopheleasonality - Seasonality and Dynamics of Indian Anophelines 📈🦟

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/326513361.svg)](https://zenodo.org/badge/latestdoi/326513361)
<!-- badges: end -->


## Overview
This repository contains the code used to analyse the results of a systematic review exploring the seasonality of various Anopheline species endemic to India. This review was carried out in order to identify entomological surveys in which mosquito collections had been conducted monthly (or finer resolution) over a period of at least a year. 

In this work, we explore the patterns of seasonality possessed by different mosquito populations, characterise these temporal dynamics in detail and develop a Bayesian framework to cluster these time-series into "dynamical archetypes" comprising groups of populations displaying similar patterns of seasonality. For example:

![alt text](https://github.com/cwhittaker1000/anopheleseasonality/blob/main/arcehtype_example.JPG?raw=true)

This work is not yet published, but has been pre-printed and is available on medRxiv: https://www.medrxiv.org/content/10.1101/2021.01.09.21249456v1.

## Repo Contents
- [Analyses](./Analyses): Code running the analyses and generates the figures featured in the paper.
- [Datasets](./Datasets): Contains entomological (from the systematic review) and environmental data (satellite-derived) used in the analyses. Collated entomological data can specifically be found in the csv  [Processed_Catch_Data.csv](./Datasets/Systematic_Review/Processed_Catch_Data.csv).
- [Figures](./Figures): Containing .PDF and .ai versions of paper figures.
- [Functions](./Functions): Extra functions required for the analyses presented in the paper.
- [Model Files](./Model_Files): STAN model files for Negative Binomial Gaussian Process fitting and Regularised Multinomial Logistic Regression. 
- [Outputs](./Outputs): Containing .rds outputs from model fitting and time-series characterisation.
- [Paper Draft](./Paper_Draft): Containing draft Manuscript and associated Supplementary Information.

## Software Requirements
Other than the required R packages (specified in each script), running the code contained in this repository requires the probabilistic programming language STAN for model fitting. STAN is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) simulation and can be accessed via the `rstan` package, available here:
- The package **rstan** (Version 2.21.2 used here) (see: https://cran.r-project.org/web/packages/rstan/index.html)
More information and details about the software and its use via R are available here: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

## Installation Guide and Instructions for Use
The following instructions require that all the relevant `R` packages have been installed by the user and that STAN has been installed. To replicate and reproduce the analyses presented in this paper, download the [Github Release](https://github.com/cwhittaker1000/anopheleseasonality/releases/tag/v1.0.0) associated with this repository. Be sure to read the intrusctions that accompany the release to ensure all the large files get placed in the right directories. Then, run the relevant set of scripts:
- Scripts in [Analyses/2_Time_Series_GP_Fitting_and_Analyses](./Analyses/2_Time_Series_GP_Fitting_and_Analyses) fit Negative Binomial Gaussian Processes to smooth the time-series, characterise their temporal properties and cluster these into dynamical archetypes.
- Scripts in [Analyses/3_Figure_Plotting](./Analyses/3_Figure_Plotting) produce the specific plots and figures present in the paper. 

## Note
- This repository contains all of the datasets/outputs generated in the analyses carried out except for a small number which are larger than GitHub's maximum allowed file size. These particularly large outputs (alongside the rest of the repo) are available for download via the pinned [Github Release](https://github.com/cwhittaker1000/anopheleseasonality/releases/tag/v1.0.0) associated with this repository. 
- The fact we have included all of these files means the repository/release is **very** large - if users do not intend to reproduce all of the analyses (but are instead interested in a particular analysis e.g. the predictive mapping or the time-series clustering), only downloading the data and files relevant to the specific analysis of interest should reduce the size significantly. 
- Scripts in [Analyses/1_Covariate_Extraction_and_Collation](./Analyses/1_Covariate_Extraction_and_Collation) carry out the raw processing of rainfall data (from CHIRPS via Google Earth Engine: https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY) and a suite of environmental of environmental covariates (sources for each detailed in [0_Raster_Processing.R](./Analyses/1_Covariate_Extraction_and_Collation/0_Raster_Processing.R)) specific to the location each study was carried out in. The total size of these raw data files is >20Gb and they are not provided with the repo - if the raw files are required, they must be redownloaded from the relevant sources. 
  - Instead, we provide the processed versions of the rainfall data (available in [Location_Specific_CHIRPS_Rainfall](./Datasets/CHIRPS_Rainfall_Data/Location_Specific_CHIRPS_Rainfall)) and environmental variables (available via the [GitHub Release](https://github.com/cwhittaker1000/anopheleseasonality/releases/tag/v1.0.0) associated with this repository) specifically required to reproduce the analyses presented here. 


Any issues, please post an issue on this Github repository or feel free to reach out at charles.whittaker16@imperial.ac.uk 😄 
