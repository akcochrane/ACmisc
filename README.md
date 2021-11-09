<!-- README.md is generated from README.Rmd. Please edit that file -->

# ACmisc

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview to ACmisc

Various functions used by Aaron Cochrane, and encouraged for
collaborators to use. Please don’t hesitate to branch, submit pull
requests, etc., or to contact Aaron directly.

## Installing the package

The R package `devtools` includes a very easy way to install packages
from Github.

    devtools::install_github('akcochrane/ACmisc', build_vignettes = TRUE)

## A couple example functions

### `BICBF()` – Bayes Factors for common models

Calculates BIC for the entire lm(), rlm(), glm(), lmer(), or glmer()
model, and for the entire set of models dropping each fixed effect one
at a time. Uses the change in BIC between models to approximate Bayes
Factors (essentially penalized Likelihood Ratios, or evidence for one
model over another). Follows the recommendations of Wagenmakers (2007;
PBR)

``` r
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 4.1.1

    ## Loading required package: Matrix

``` r
mod <- lmer(mot ~ nback * ravens + (1|isAdult),dat_cochraneEtAl_2019_PLOSOne)

BF_table <- BICBF(mod)
```

    ## Loading required package: MASS

``` r
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
the density estimation function from the RWiener package. First fits a
standard 4-parameter Wiener model to the entire RT+response-category
vector. Then uses the boundary separation, bias, and non-decision time
from this overall model, and finds the best set of parameters to create
a drift rate vector as a linear

### `findLowerRT()`

Follows the procedure of Ratcliff & Tuerlinckx (2002) to find the lowest
RT at which accuracy is below chance.
