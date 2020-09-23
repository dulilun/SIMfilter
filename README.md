
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SIMfilter

<!-- badges: start -->

<!-- badges: end -->

The goal of SIMfilter is to provide an R package for our SIM paper; see
the link for detials: <https://projecteuclid.org/euclid.aos/1403715201>

## Installation

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dulilun/SIMfilter")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SIMfilter)
## basic example code
p_1 = rnorm(100)
p_1[1:10]=p_1[1:10]+2
p_1 = 1-pnorm(p_1)
p_2 = rnorm(100)
p_2[1:10]=p_2[1:10]+2
p_2 = 1-pnorm(p_2)
out = SIM(p_1, p_2, alpha = 0.1, SIM_F_0_method=1, method_theta=1, method_lambda=1, method_t=3)
```
