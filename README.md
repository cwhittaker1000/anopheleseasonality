# anopheleasonality - Seasonality and Dynamics of Indian Anophelines 📈🦟

## Overview
This repository contains the code used to analyse the results of a systematic review exploring the seasonality of various Anopheline species endemic to India. This review was carried out in order to identify entomological surveys in which mosquito collections had been conducted monthly (or finer resolution) over a period of at least a year. 

In this work, we explore the patterns of seasonality possessed by different mosquito populations, characterise these temporal dynamics in detail and develop a Bayesian framework to cluster these time-series into "dynamical archetypes" comprising groups of populations displaying similar patterns of seasonality. This work is not yet published, but has been pre-printed and is available on medRxiv: https://www.medrxiv.org/content/10.1101/2021.01.09.21249456v1.

## Repo Contents
- [Analyses](./Analyses): Code running the analyses and generates the figures featured in the paper.
- [Datasets](./Datasets): Contains entomological (from the systematic review) and environmental data (satellite-derived) used in the analyses.
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
The following instructions require that all the relevant `R` packages have been installed by the user and that STAN has been installed. To replicate and reproduce the analyses presented in this paper, clone this repository and download it your local machine. Then, run the relevant set of scripts:
- Scripts in [Analyses/1_Covariate_Extraction_and_Collation](./Analyses/1_Covariate_Extraction_and_Collation) extract rainfall data and a suite of environmental covariates specific to the location each study was carried out in. 
- Scripts in [Analyses/2_Time_Series_GP_Fitting_and_Analyses](./Analyses/2_Time_Series_GP_Fitting_and_Analyses) fit Negative Binomial Gaussian Processes to smooth the time-series, characterise their temporal properties and cluster these into dynamical archetypes.
- Scripts in [Analyses/3_Figure_Plotting](./Analyses/3_Figure_Plotting) produce the specific plots and figures present in the paper. 

Any issues, please post an issue on this Github repository or feel free to reach out at charles.whittaker16@imperial.ac.uk 😄
