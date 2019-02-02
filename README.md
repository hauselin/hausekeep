
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

An example of the `es()` effect size conversion function.

``` r
library(hausekeep) # load package
es(d = 0.2)
#> d: 0.2
#>     d   r   R2   f oddsratio logoddsratio   auc
#> 1 0.2 0.1 0.01 0.1     1.437        0.363 0.579
```

Example plot
------------

<img src="man/figurespressure-1.png" width="100%" />
