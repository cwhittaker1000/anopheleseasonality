# Indian Anopheline Seasonality Systematic Review

## Overview
This repository contains the code used to analyse the results of a systematic review exploring the seasonality of various Anopheline species endemic to the Indian subcontinent. Briefly, a systematic review was carried out in order to identify entomological surveys in which mosquito collections had been conducted monthly (or finer resolution) over a period of at least a year. The collated entomological data was then analysed using a Bayesian Gaussian Process based approach in order to explore the patterns of seasonality displayed by different mosquito species.

## Repo Contents
- [Analyses](./Analyses): Code running the analyses and generates the figures featured in the paper.
- [Functions](./Functions): Extra functions required for the analyses presented in the paper.
- [Data](./Data): Containing the individual-level social contact data collated as part of the systematic review and that are used in the analyses presented here. 
- [Figures](./Figures): Containing .PDF and .ai versions of paper figures.
- [Outputs](./Outputs): Containing .rds outputs from model fitting (undertaken using a Bayesian Regression approach and implemented  using BRMS) for each of the analyses conducted.

## Software Requirements
Running the code contained in this repository requires the following
- The R Programming Language (Version 4.0.2 "Taking Off Again" used here) 
- The package **rstan** (Version 2.21.2 used here) (see: https://cran.r-project.org/web/packages/rstan/index.html)
- The package **brms** (Version 2.13.5 used here) (see: https://cran.r-project.org/web/packages/brms/index.html)
- The packages contained within **tidyverse** (Version 2.13.5 used here) (see: https://cran.r-project.org/web/packages/tidyverse/index.html)
- The package **gdata** (version 2.13.5 used here) (see: https://cran.r-project.org/web/packages/gdata/index.html)

This work utilises the probabilistic programming language STAN for model fitting (implemented via the package BRMS). STAN is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) simulation. More information and details about the software and its use via R are available here: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

## Installation Guide and Instructions for Use
The following instructions require that all the relevant `R` packages have been installed by the user and that STAN has been installed. To replicate and reproduce the analyses presented in this paper, do the following: 

1. Download the [Data](./Data) folder of this repository. 
2. Download the `R` code from  [Analyses](./Analyses) for the particular part of the analysis you are trying to reproduce. 
4. Run the `R` code. The output from running this code will be a number of MCMC objects, as well as a series of plots representing the output from MCMC based fitting of the relevant model to the collated data. These plots form the basis of the figures presented in the publication. 
