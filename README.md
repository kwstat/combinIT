

<!-- README.md is generated from README.Rmd. Please edit that file -->

# combinIT
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

There are several non-functional-form-based interaction tests for testing interaction in unreplicated two-way layouts. However, no single test can detect all patterns of possible interaction and the tests are sensitive to a particular pattern of interaction. This package combines six non-functional-form-based interaction tests for testing additivity. These six tests were proposed by Boik (1993a), Piepho (1994), Kharrati-Kopaei and Sadooghi-Alvandi (2007), Franck et al. (2013), Malik et al. (2016) and Kharrati-Kopaei and Miller (2016). The p-values of these six tests are combined by Bonferroni, Sidak, Jacobi polynomial expansion, and the Gaussian copula methods to provide researchers with a testing approach which leverages many existing methods to detect disparate forms of non-additivity. This package is based on the following published paper: Zahra Shenavari & Mahmood Kharrati‐Kopaei, 2018. "A Method for Testing Additivity in Unreplicated Two‐Way Layouts Based on Combining Multiple Interaction Tests," International Statistical Review, International Statistical Institute, vol. 86(3), pages 469-487, December. In addition, several sentences in help files or descriptions were copied from that paper.

# README Notes

The reader should note that the p-value of interaction tests is obtained by a Monte Carlo simulation procedure in most cases. In addition, the applied tests are based on the normality assumption of data. To calculate the p-value of a test very fast, several codes have been written in Rcpp. However, the speed of running the codes depends on the number of Monte Carlo samples, nsim. Note that the most Monte Carlo error of estimating a p-value is 1.96*sqrt(1/(4*nsim)) with probability 0.95. Finally, the reader should note that this package has multiple package dependencies. 


# Installation

You can install combinIT from github with:

``` r
# install.packages("devtools")
devtools::install_github("haghbinh/combinIT")
```

