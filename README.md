# KMSim: Simulating Survival Times from a Kaplan-Meier Estimator (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`KMSim` is an R package for simulating times to event from a
**Kaplan-Meier (KM) estimator** of a survival function, constructed from an
observed sample of (possibly censored) survival times.

The standard KM estimator defines a discrete, non-parametric estimate of the
survival function. `KMSim` uses this estimated distribution to draw new
survival times via inverse transform sampling, enabling non-parametric
simulation without fitting a parametric model. Simulation can be performed
for the entire sample, a subgroup defined by covariates, or a single
individual, making the package useful for bootstrapping, sample size
calculation, and simulation studies in survival analysis.

## Installation

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/KMSim")
library(KMSim)
```

## Main function

| Function | Description |
|---|---|
| `KMSim` | Simulate survival times from a Kaplan-Meier estimator |

The function takes a `survfit` object (from the `survival` package) as input
and returns simulated times to event. For full documentation: `?KMSim`

## Quick example

```r
library(KMSim)
library(survival)

# Fit a KM estimator to the lung dataset
km_fit <- survfit(Surv(time, status) ~ 1, data = lung)

# Simulate 100 new survival times from the estimated distribution
sim_times <- KMSim(km_fit, n = 100)
```

## Tutorial

- [KMSim: Simulating from a Kaplan-Meier estimator](https://rpubs.com/FJRubio/KMSim) —
  illustrative examples using real survival data

## Related resources

- [SimLT](https://github.com/FJRubio67/SimLT) — simulation of survival times
  from a life table (parametric piecewise-constant hazard approach)
- [HazReg](https://github.com/FJRubio67/HazReg) — parametric hazard-based
  regression models for survival data
- [survival](https://cran.r-project.org/web/packages/survival/index.html) —
  the R package used to produce the `survfit` input object

## License

This package is licensed under the [MIT License](LICENSE).
