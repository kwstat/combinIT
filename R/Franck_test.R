#' Franck's (2013) et al. Test for Interaction
#'
#' This function calculates Franck's (2013) et al. test statistic, ACMIF, and corresponding p-value.
#'
#' @param x numeric matrix, \eqn{a \times b} data matrix where the number of row and column is corresponding to the number of factor levels.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param Elapsed.time logical: if \code{TRUE} the progress will be printed in the console.
#' @param alpha a numeric value, the level of the test. The default value is 0.05.
#'
#' @details Franck et al. (2013) derived a test statistic based on the “hidden additivity” structure.
#'  They defined this structure as “the levels of one factor belong in two or more groups such that within each group the effects of the two factors are additive but the groups may interact with the ungrouped factor”.
#'  To detect hidden additivity, Franck et al. (2013) divided the table of data into two sub-tables (based on the rows of the data matrix) and an interaction F-test was developed.
#'  Then, they performed a search over all possible configures of data and used the maximum of the interaction F-test as a test statistic. The hypothesis of no interaction is rejected when the maximum interaction F-test is large.
#'  Note that the number of rows should be greater than two. This test is powerful when there is a hidden additivity structure in the data set.
#'
#'
#' @return An object of the class \code{ITtest}, which is a list inducing following components::
#' \item{pvalue.exact}{The calculated exact Monte Carlo p-value.}
#' \item{pvalue.appro}{The Bonferroni-adjusted p-value is calculated.}
#' \item{statistic}{The value of the test statistic.}
#' \item{Nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#' \item{Level}{The level of test.}
#' \item{Result}{The result of the test at the alpha level with some descriptions on the type of significant interaction.}
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
#' data(CNV)
#' Franck.test(CNV, nsim = 1000, Elapsed.time = FALSE)
#'
#' @importFrom stats pchisq pf qnorm var
#' @export
Franck.test <- function(x, nsim = 10000, alpha = 0.05, Elapsed.time = TRUE) {
  DNAME <- deparse1(substitute(x))
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    bl <- nrow(x)
    tr <- ncol(x)
    n <- tr * bl
    if (bl < 3) {
      warning("Frank.test needs at least three levels for the row factor.")
      str <- Result.Franck(x, nsim = nsim, alpha = alpha, simu = NULL)
      out <- list(
        pvalue.exact = NA,
        pvalue.appro = NA,
        nsim = nsim,
        statistic = NA,
        data.name = DNAME,
        test = "Franck Test",
        Level = alpha,
        Result = str
      )
    } else {
      cch <- 2^(bl - 1) - 1
      statistics <- hh_f(x)
      simu <- rep(0, 0)
      if (Elapsed.time) {
        pb <- completed(nsim)
        for (i in 1:nsim) {
          simu[i] <- hh_f(matrix(rnorm(n), nrow = bl, ncol = tr))
          if (i == pb$pr[pb$j]) pb <- nextc(pb, i)
        }
      } else {
        for (i in 1:nsim) {
          simu[i] <- hh_f(matrix(rnorm(n), nrow = bl, ncol = tr))
        }
      }
      hidden <- mean(statistics < simu)
      adjpvalue <- (1 - pf(statistics, (tr - 1), (tr - 1) * (bl - 2))) * cch
      hidden.apr <- min(1, adjpvalue)
      qFranck <- quantile(simu, prob = 1 - alpha, names = FALSE)
      if (hidden < alpha) {
        str <- Result.Franck(x, nsim = nsim, alpha = alpha, simu = simu)
      } else {
        str <- paste("The Franck.test could not detect any significant interaction.", "The estimated critical value of the Franck.test with", nsim, "Monte Carlo samples is:", round(qFranck, 4), "\n")
      }
      out <- list(
        pvalue.exact = hidden,
        pvalue.appro = hidden.apr,
        nsim = nsim,
        statistic = statistics,
        data.name = DNAME,
        test = "Franck Test",
        Level = alpha,
        Result = str
      )
    }
    structure(out, class = "ITtest")
  }
}
