<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/425863127.svg)](https://zenodo.org/badge/latestdoi/425863127)

# ACmisc

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview to ACmisc

Various functions used by Aaron Cochrane, and encouraged for
collaborators to use. Please don’t hesitate to branch, submit pull
requests, etc., or to contact Aaron directly. And if you use the
functions in published work, it would be wonderful if you cited this
package using the Zenodo DOI. Thanks!

## Installing the package

The R package `devtools` includes a very easy way to install packages
from Github.

    devtools::install_github('akcochrane/ACmisc', build_vignettes = TRUE)

## A handful of example functions

### `BICBF()` – Bayes Factors for common models

Calculates BIC for an entire `lm()`, `rlm()`, `glm()`, `lmer()`, or
`glmer()` model, and for the entire set of models dropping each fixed
effect one at a time. Uses the change in BIC between models to
approximate Bayes Factors (essentially penalized Likelihood Ratios, or
evidence for one model over another). Follows the recommendations of
Wagenmakers (2007; PBR).

``` r
library(lme4)
mod <- lmer(mot ~ nback * ravens + (1|isAdult),dat_cochraneEtAl_2019_PLOSOne)

BF_table <- BICBF(mod)
knitr::kable(BF_table,digits=3)
```

|              | Estimate | Std..Error | t.value | BIC_dropped | BFlog3 |
|:-------------|---------:|-----------:|--------:|------------:|-------:|
| (Intercept)  |    0.904 |      0.186 |   4.864 |    -316.390 |     NA |
| nback        |   -0.105 |      0.221 |  -0.475 |    -316.859 | -0.213 |
| ravens       |   -1.084 |      0.412 |  -2.632 |    -310.396 |  2.728 |
| nback:ravens |    1.289 |      0.479 |   2.691 |    -310.054 |  2.884 |

### `d1iso()`

Extracts scores from a single isotonic \[think: rank-based\] underlying
dimension from multivariate data. Very useful for extracting the common
variation in data when other methods (e.g., PCA) would struggle due to
different variables’ lack of Gaussian-ness.

### `ddm_dr_lm()` – Drift Diffusion Linear Model

Essentially fits a generalized linear model for a Wiener process, using
the density estimation function from the `RWiener` package. First fits a
standard 4-parameter Wiener model to an entire vector of RTs and
response-categories. Then uses the boundary separation, bias, and
non-decision time from this overall model, and finds the best set of
parameters to create a drift rate vector as a linear function of the
right-hand side of the model formula. **Only numeric predictors have
been tested.**

### `findLowerRT()`

Follows the procedure of Ratcliff & Tuerlinckx (2002) to find the lowest
RT at which accuracy is below chance.

### `getTime()`

Returns a string with a formatted date and time, to the nearest minute.

### `isOutlier()`

Uses robust covariance estimation of Mahalanobis distances to identify
multivariate outliers within a dataset; follows the procedure outlined
by Leys and colleagues (2018; JESP).

### `resetSeed()`

Re-“randomizes” the seed for random generation.

### `robustLM_bayes()`

Implements a very robust linear model, with mixed-effects terms
possible. Uses median regression within the `brms` package to fit the
specified model and sub-models, and uses the `bridgesampling` package to
estimate Bayes Factors supporting the model *with* the predictor
compared to the model *without* that predictor.

### `sourceIfChanged()`

Useful when scripts may take a long time to run. Runs a file, just like
`source()`. But the resulting namespace is saved (in an `.RData`), as
well as a plain-text copy. The next time that `sourceIfChanged` is run,
it simply loads the `.RData` file if there haven’t been any changes, and
it re-`source`s the `.R` file if there has been changes.

### `YeoJohn()`

Applies a Yeo-Johnson univariate power transformation to a variable,
after optimizing the Yeo-Johnson parameter to minimize the variable’s
skew.
