#' Kharrati-Kopaei and Sadooghi-Alvandi's (2007) test for interaction
#'
#' This function calculates Kharrati-Kopaei and Sadooghi-Alvandi's test statistic and corresponding p-value for testing interaction.
#'
#' @param x numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
#'   , respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param Elapsed.time logical: if \code{TRUE} the progress will be printed in the console.
#' @details  Suppose that \eqn{a \le 3} and \eqn{a \ge 4}. Consider the \eqn{l}-th division of the data table into two sub-tables,
#'  obtained by putting \eqn{b_1} (\eqn{2<U+2264>b_1<U+2264>b-2}) rows in the first sub-table and the remaining \eqn{b_2} rows in the second sub-table (\eqn{b_1+b_2=a}).
#'  Let RSS1 and RSS2 denote the residual sum of squares for these two sub-tables, respectively. For a particular division \eqn{l}, let \eqn{F_l=max<U+2061>(F_l,1/F_l }
#'  where \eqn{F_l=(b_2-1)RSS1/((b_1-1)RSS2)} and let \eqn{P_l} denote the corresponding p-value.
#'  Kharrati-Kopaei and Sadooghi-Alvandi (2007) proposed their test statistic as the minimum value of \eqn{P_l} over \eqn{l=1,…,2^(b-1)-b-1} all possible divisions of the table.
#'  Note that if the rows number, \eqn{b}, of data matrix is less than the columns number, \eqn{a}, the data matrix is transposed. In addition, this method of testing requires that the data matrix has more than three
#'  rows or columns. This test procedure is powerful for detecting interaction when the magnitude of interaction effects is heteroscedastic across the sub-tables of observations.

#' @return An object of the class \code{ITtest}, which is a list inducing following components::
#' \item{pvalue.exact}{The calculated exact Monte Carlo p-value.}
#' \item{pvalue.appro}{The Bonferroni-adjusted p-value is calculated.}
#' \item{statistic}{The value of the test statistic.}
#' \item{Nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#'
#'
#' @references Kharrati-Kopaei, M., Sadooghi-Alvandi, S.M. (2007). A New Method for
#'  Testing Interaction in Unreplicated Two-Way Analysis of Variance. Communications
#'  in Statistics-Theory and Methods 36:2787–2803.
#'
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#'
#' @examples
#' data(IDCP)
#' KKSA.test(IDCP, nsim = 1000, Elapsed.time = FALSE)
#' 
#' @export
KKSA.test <- function(x, nsim = 10000, Elapsed.time = TRUE) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    DNAME <- deparse1(substitute(x))
    bl <- nrow(x)
    tr <- ncol(x)
    n <- tr * bl
    if (bl < tr) {
      warning("The input data matrix is trasposed")
      x <- t(x)
      te <- bl
      bl <- tr
      tr <- te
    }
    if (bl < 4) {
      warning("KKSA needs at least 4 levels for a factor")
      out <- list(
        pvalue.exact = NA,
        pvalue.appro = NA,
        nsim = nsim,
        statistic = NA,
        data.name = DNAME,
        test = "KKSA Test")
    } else {
      cck <- 2^(bl - 1) - 1 - bl
      statistics <- kk_f(x)
      simu <- rep(0, 0)
      if (Elapsed.time) {
        pb <- completed(nsim)
        for (i in 1:nsim) {
          simu[i] <- kk_f(matrix(rnorm(n), nrow = bl))
          if (i == pb$pr[pb$j]) pb <- nextc(pb, i)
        }
      } else {
        for (i in 1:nsim) {
          simu[i] <- kk_f(matrix(rnorm(n), nrow = bl))
        }
      }
      KKSA.p <- mean(statistics > simu)
      KKSA.p.apr <- statistics * cck
      KKSA.p.apr <- min(1, KKSA.p.apr)
      out <- list(
        pvalue.exact = KKSA.p,
        pvalue.appro = KKSA.p.apr,
        nsim = nsim,
        statistic = statistics,
        data.name = DNAME,
        test = "KKSA Test"
      )
    }
    structure(out, class = "ITtest")
  }
}
