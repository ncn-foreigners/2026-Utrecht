# 2026-Utrecht

Simulation study comparing non-probability survey estimators using quantile-based covariate representations. Prepared for the 2026 Utrecht conference.

## Overview

The study evaluates IPW, mass imputation (MI), and doubly robust (DR) estimators from the [`nonprobsvy`](https://github.com/ncn-foreigners/nonprobsvy) package under different covariate basis specifications:

- **main** -- original covariates (x1, x2)
- **quartiles** -- cumulative quantile indicators at quartiles
- **main+quartiles** -- original covariates + quartile indicators
- **deciles** -- cumulative quantile indicators at deciles
- **main+deciles** -- original covariates + decile indicators

Two non-probability samples with different selection mechanisms (bd1, bd2) are drawn from a simulated population (N=20,000) with both continuous (y11, y12) and binary (y21, y22) outcomes. Each simulation runs 500 replications using parallel `foreach` with `doRNG` for reproducibility.

### Simulation study 1 (`sim-study.R`)

Standard polynomial DGP with linear (y11) and quadratic (y12) continuous outcomes, plus corresponding binary outcomes (y21, y22). Selection mechanisms: bd1 uses monotone selection via x2, bd2 uses non-linear selection via squared deviations of x1 and x2.

### Simulation study 2 (`sim-study-2.R`)

Smooth non-monotone DGP designed to demonstrate the advantage of quantile-based representations. Outcomes are continuous functions of the empirical CDF of x2 (Fx2 = rank(x2)/N):

- **y11**: U-shaped -- `4(2Fx2 - 1)^2 + 0.3x1` (high at tails, low at median)
- **y12**: varying slopes -- `1 + (3Fx2 - 1.5) * x1` (slope of x1 changes sign across the distribution of x2)
- **y21, y22**: binary analogues via logistic link (U-shaped and inverted-U probabilities)

Selection: bd1 oversamples the middle of x2 (IQR region), bd2 oversamples the extremes. Smooth functions ensure that quantile estimation error from the finite probability sample does not cause boundary-related bias, and finer quantile grids (deciles) progressively improve the piecewise-constant approximation.

## Repository structure

```
codes/
  functions.R        -- helper functions (make_quantile_basis, extract, check_bin_counts)
  sim-study.R        -- simulation study 1 (polynomial DGP)
  sim-study-2.R      -- simulation study 2 (smooth non-monotone DGP)
  reporting.qmd      -- Quarto report for study 1 (bias, SE, RMSE, CI coverage)
  reporting-2.qmd    -- Quarto report for study 2
results/
  sim-results.rds    -- simulation output from study 1
  sim-results-2.rds  -- simulation output from study 2
```

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/) to pin exact package versions. To restore the environment:

```r
# install renv if needed
install.packages("renv")

# restore all packages from renv.lock
renv::restore()
```

## How to run

```r
# run simulation study 1 (saves results/sim-results.rds)
Rscript codes/sim-study.R

# run simulation study 2 (saves results/sim-results-2.rds)
Rscript codes/sim-study-2.R

# render reports (produces self-contained HTML + PNGs in results/figures/)
quarto render codes/reporting.qmd
quarto render codes/reporting-2.qmd
```

## R packages

survey, data.table, nonprobsvy, doSNOW, progress, foreach, doRNG, ggplot2
