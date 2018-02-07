# Simulations of GIV regression 

## Introduction
This repository contains the scripts for the SNP level simulations from [DiPrete et al. (2018) *Genetic Instrumental Variable (GIV) regression: Explaining socioeconomic and health outcomes in non-experimental data*](https://www.biorxiv.org/content/early/2018/02/05/134197.article-info). See that text and the supporting information for the technical details on these simulations.

All scripts have been written in R. 

## Overview
This repository contains six R scripts, each script tests various methods (including OLS, MR and GIV) in a different setting. 

#### heritability.r
This simulation shows how polygenic scores can be used to estimate narrow-sense heritability. See SI section 1 and SI table 1.

#### pleiotropy.r
In this script we simulate a model where the goal is to estimate the effect of an exposure on an outcome. This model has endogeneity due to pleiotropy between the exposure and the outcome. See SI section 2 and SI tables 2-6.

#### pleiotropy_table7.r
This script is designed specifically for SI table 7. The parameter settings were unobtainable in `pleiotropy.r`.

#### genetic_endogeneity.r
Here we simulate a model where next to pleiotropy there is additional endogeneity that is genetic-related (e.g. genetic nurturing). See SI section 3 and SI tables 8-10.

#### unrelated_endogeneity.r
In this script a model is simulated where the additional endogeneity is not related to any genetic factors. See SI section 4 and SI tables 11-12.

#### partial_control.r
In this script a model is simulated where the additional endogeneity is not related to any genetic factors, but we can partially control for this endogeneity in the hold-out sample. See SI section 4 and SI tables 13-15.
