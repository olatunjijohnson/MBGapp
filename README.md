
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MBGapp

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/olatunjijohnson/MBGapp.svg?branch=main)](https://travis-ci.org/olatunjijohnson/MBGapp)
[![Codecov test
coverage](https://codecov.io/gh/olatunjijohnson/MBGapp/branch/master/graph/badge.svg)](https://codecov.io/gh/olatunjijohnson/MBGapp?branch=main)
<!-- badges: end -->

The goal of MBGapp is to allow user to explore and analyse
geostatistical
data

## Installation

<!-- You can install the released version of MBGapp from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("MBGapp") -->

<!-- ``` -->

You can install the development version on github with

``` r
devtools::install_github("olatunjijohnson/MBGapp", ref="main")
```

## Example

Download an example dataset from
[here](https://drive.google.com/uc?export=download&id=14MPEAqI7qIP-U9q_vbGuDG8hkoFSjc9A).
The data is the Loaloa prevalence data in Cameroon.

This is a basic example which shows you how to run the app:

``` r
library(MBGapp)
## run the App
run_app()  # use the code
```

## Alternative way to run in R

You can also run the following line of code to run in
R

``` r
shiny::runGitHub(repo="MBGapp", username= "olatunjijohnson", ref="main", subdir = "inst/MBGapp")
```

## Online version

The app can also be accessed online via the following link:

<https://olatunjijohnson.shinyapps.io/mbgapp/>
