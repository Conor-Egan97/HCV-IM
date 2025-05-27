# HCV-IM: Modelling Hepatitis C Incidence in People Who Inject Drugs

This repository contains R scripts and documentation for modelling the incidence of Hepatitis C Virus (HCV) among people who inject drugs (PWID). It uses simulated data representative of the Unlinked Anonymous Monitoring (UAM) Survey and mirrors the methodology detailed in the publication _"Analysing HCV incidence trends in people who inject drugs using serial behavioural and seroprevalence data: A modelling study"_.

## Project Overview

The aim of this project is to model HCV incidence trends in PWID using behavioural and serological survey data collected over time. Three statistical models are applied to the data to estimate the force of infection (FOI) across calendar time and injecting duration. Confidence intervals are generated using bootstrap methods.

This project facilitates reproducibility of the modelling approaches used in the paper and supports further public health interpretation and visualisation of HCV incidence dynamics.

## Repository Structure

- `sim_UAM.R`: Simulates synthetic UAM-like data that mimics the structure of the real UAM dataset. This allows public users to explore and test the models in the absence of access to real survey data.
- `PC_FOI.R`: Fits a multiplicative **piecewise constant** model without interactions between injecting duration and calendar time. Confidence intervals for the FOI are derived using bootstrap methods.
- `FP_FOI.R`: Fits a **fractional polynomial** model with interaction terms between injecting duration and calendar time. The model is applied to the cumulative hazard to derive the FOI.
- `NCS_FOI.R`: Fits a **natural cubic spline** model to the cumulative hazard, incorporating spline interactions in both calendar time and injecting duration dimensions. Confidence intervals are also derived via bootstrapping.
- `PC_functions/`: Contains four R scripts required to support the piecewise constant model:
  - `fn_getparamATnoint.R`
  - `fn_mkcut.R`
  - `fn_pred.R`
  - `fn_univarATnoint.R`  
  These scripts were written by Dr. Ross Harris, Senior Statistician at the UK Health Security Agency (UKHSA).

## Data

The real UAM data are not publicly available. This repository instead uses a simulated version created via `sim_UAM.R` to mimic key characteristics of the real dataset.

For those interested in the UAM survey, the UKHSA publishes an annual report with relevant details.  
[Unlinked Anonymous Monitoring Survey of HIV and Viral Hepatitis among PWID – 2024 Report](https://www.gov.uk/government/publications/people-who-inject-drugs-hiv-and-viral-hepatitis-monitoring/unlinked-anonymous-monitoring-uam-survey-of-hiv-and-viral-hepatitis-among-people-who-inject-drugs-pwid-2024-report)

## Modelling Approaches

### 1. Piecewise Constant (PC) Model
- Multiplicative model with stratified calendar time and injecting duration intervals.
- No interaction terms included.
- Bootstrap-derived confidence intervals.
- Requires external support functions (`PC_functions/`).

### 2. Fractional Polynomial (FP) Model
- Fits second-degree fractional polynomials of injecting duration.
- Includes interactions between injecting duration and calendar time.
- Model fitted to cumulative hazard; FOI is derived via differentiation.
- Bootstrap used to construct FOI confidence intervals.

### 3. Natural Cubic Spline (NCS) Model
- Applies 7 knots (injecting duration) and 3 knots (calendar time).
- Models cumulative hazard using spline basis functions with interaction terms.
- FOI computed from the fitted hazard model.
- Confidence intervals estimated via bootstrapping.

## Outputs

Each script produces:
- Estimates of the force of infection (FOI) by injecting duration and calendar time.
- Bootstrap confidence intervals around the FOI estimates.
- ggplot2-based visualisations of FOI trends and modelled uncertainty.

## Requirements

This project uses the following R packages:

```r
install.packages(c(
  "magrittr", "dplyr", "ggplot2", "boot", "bbmle", "stats4"
))
