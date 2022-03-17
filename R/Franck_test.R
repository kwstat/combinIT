#' Franck's (2013) et al. test for interaction
#'
#' This function calculates Franck (2013) et al. test statistic,ACMIF, and corresponding p-value.
#'
#' @param x numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
#'   , respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param dist a character, if dist="sim", a Monte Carlo simulation is used for calculating exact p-value,
#'  and if dist="adj", the Bonferroni-adjusted p-value is calculated. The default is "sim".
#'
#' @details Franck et al. (2013) derived a test statistic based on the “hidden additivity” structure.
#'  They defined this structure as “the levels of one factor belong in two or more groups such that within each group the effects of the two factors are additive but the groups may interact with the ungrouped factor”.
#'  To detect hidden additivity, Franck et al. (2013) divided the table of data into two sub-tables and an interaction F-test was developed.
#'  Then, they performed a search over all possible configures of data and used the maximum of the interaction F-test as a test statistic. The hypothesis of no interaction is rejected when the maximum interaction F-test is large.
#'  Note that, if rows number, \eqn{b}, of data matrix is less than the columns number, \eqn{a},
#'  the data matrix is transposed. Note that the this test method is powerful when there is a hidden additivity structure in the data set.
#'
#' @return pvalue, calculated p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#'
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#'
#' @references
#'  Franck, C., Nielsen, D., Osborne, J.A. (2013). A method for detecting hidden additivity in two-factor unreplicated experiments.
#'  Computational Statistics and Data Analysis 67:95-104.
#'
#'  Franck, C., Osborne, J.A. (2016).  Exploring Interaction Effects in Two-Factor Studies using the hidden Package in R.
#'  R Journal 8 (1):159-172.
#'
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#'
#' @examples
#' \dontrun{
#' data(CNV)
#' Franck.test(CNV, nsim = 10000, dist = "sim")
#' }
#' @importFrom stats pchisq pf qnorm var
#' @export
Franck.test <- function(x, nsim = 10000, dist = "sim") {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    bl <- nrow(x)
    tr <- ncol(x)
    n <- tr * bl
    if (bl < tr) {
      warning("The input data matrix is transposed")
      x <- t(x)
      te <- bl
      bl <- tr
      tr <- te
    }
    cch <- 2^(bl - 1) - 1
    statistics <- hh_f(x)
    if (dist != "sim" & dist != "adj") stop("\"dist\" parameter should be equal to \"sim\" or \"adj\".")
    if (dist == "sim") {
      simu <- rep(0, 0)
      for (i in 1:nsim) {
        simu[i] <- hh_f(matrix(rnorm(n), nrow = bl, ncol = tr))
        cat(paste(round(i / nsim * 100), "% completed"), "\n")
        if (i == nsim) {
          cat(": Done", "\n")
        } else {
          cat("\014", "\n")
        }
      }
      hidden <- mean(statistics < simu)
    }
    if (dist == "adj") {
      adjpvalue <- (1 - pf(statistics, (tr - 1), (tr - 1) * (bl - 2))) * cch
      hidden <- min(1, adjpvalue)
    }
    out <- list(
      pvalue = hidden,
      nsim = nsim,
      dist = dist,
      statistic = statistics
    )
    return(out)
  }
}


