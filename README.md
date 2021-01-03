# Indian Anopheline Seasonality Systematic Review

## Overview
This repository contains the code used to analyse the results of a systematic review exploring the seasonality of various Anopheline species endemic to the Indian subcontinent. Briefly, a systematic review was carried out in order to identify entomological surveys in which mosquito collections had been conducted monthly (or finer resolution) over a period of at least a year. The collated entomological data was then analysed using a Bayesian Gaussian Process based approach in order to explore the patterns of seasonality displayed by different mosquito species.

## Repo Contents
- [Analyses](./Analyses): Code running the analyses and generates the figures featured in the paper.
- [Conference Presentations](./Conference Presentations): Containing presentations in preparation for various conferences. 
- [Figures](./Figures): Containing .PDF and .ai versions of paper figures.
- [Functions](./Functions): Extra functions required for the analyses presented in the paper.
- [Model Files](./Model Files): STAN model files for Negative Binomial Gaussian Process fitting and Penalised Multinomial Logistic Regression. 
- [Outputs](./Outputs): Containing .rds outputs from model fitting and time-series characterisation.
- [Paper Draft](./Paper_Draft): Containing draft Manuscript and associated Supplementary Information.

## Software Requirements
Other than the required R packages, running the code contained in this repository requires the  probabilistic programming language STAN for model fitting. STAN is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) simulation and can be accessed via the `rstan` package, available here:
- The package **rstan** (Version 2.21.2 used here) (see: https://cran.r-project.org/web/packages/rstan/index.html)
More information and details about the software and its use via R are available here: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

## Installation Guide and Instructions for Use
The following instructions require that all the relevant `R` packages have been installed by the user and that STAN has been installed. To replicate and reproduce the analyses presented in this paper, do the following: 

1. Download the [Data](./Data) folder of this repository. 
2. Download the `R` code from  [Analyses](./Analyses) for the particular part of the analysis you are trying to reproduce. 
4. Run the `R` code. The output from running this code will be a number of MCMC objects, as well as a series of plots representing the output from MCMC based fitting of the relevant model to the collated data. These plots form the basis of the figures presented in the publication. 
