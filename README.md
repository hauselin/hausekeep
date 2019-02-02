
<!-- README.md is generated from README.Rmd. Please edit that file -->
hausekeep
=========

[![Travis build status](https://travis-ci.org/hauselin/hausekeep.svg?branch=master)](https://travis-ci.org/hauselin/hausekeep) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/hauselin/hausekeep?branch=master&svg=true)](https://ci.appveyor.com/project/hauselin/hausekeep) [![Coverage status](https://codecov.io/gh/hauselin/hausekeep/branch/master/graph/badge.svg)](https://codecov.io/github/hauselin/hausekeep?branch=master)

[![DOI](https://zenodo.org/badge/168783741.svg)](https://zenodo.org/badge/latestdoi/168783741)

The goal of hausekeep is to ...

Installation
------------

To install the package, type the following commands into the R console:

``` r
# install.packages("devtools")
devtools::install_github("hauselin/hausekeep") # you might have to install devtools first (see above)
```

Examples
--------

#### `es()` function to convert between effect sizes

The `es` function converts one effect size into other effect sizes (e.g., d, r, R<sup>2</sup>, f, odds ratio, log odds ratio, area-under-curve \[AUC\]). Note that AUC calculations are slightly off!

``` r
library(hausekeep) # load package
es(d = 0.2)
#> d: 0.2
#>     d   r   R2   f oddsratio logoddsratio   auc
#> 1 0.2 0.1 0.01 0.1     1.437        0.363 0.579
es(r = c(0.1, 0.4, 0.7))
#> r: 0.1 r: 0.4 r: 0.7
#>       d   r   R2     f oddsratio logoddsratio   auc
#> 1 0.201 0.1 0.01 0.101     1.440        0.365 0.580
#> 2 0.873 0.4 0.16 0.436     4.871        1.583 0.809
#> 3 1.960 0.7 0.49 0.980    35.014        3.556 0.975
```

#### `fit_ezddm()` function to fit EZ-diffusion model for two-choice response time tasks

``` r
library(rtdists) # load package to help us simulate some data
data1 <- rdiffusion(n = 100, a = 2, v = 0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
data2 <- rdiffusion(n = 100, a = 2, v = -0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
dataAll <- rbind(data1, data2) # join data
dataAll$response <- ifelse(dataAll$response == "upper", 1, 0) # convert responses to 1 and 0
dataAll$subject <- rep(c(1, 2), each = 100) # assign subject id
dataAll$cond1 <- sample(c("a", "b"), 200, replace = T) # randomly assign conditions a/b
dataAll$cond2 <- sample(c("y", "z"), 200, replace = T) # randomly assign conditions y/z

# fit model to just entire data set (assumes all data came from 1 subject)
fit_ezddm(data = dataAll, rts = "rt", responses = "response")
# fit model to each subject (no conditions)
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject") 
# fit model to each subject by cond1
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1") 
# fit model to each subject by cond1,cond2
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2"))
```
