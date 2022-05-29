#' Piepho's (1994) Test for Interaction
#'
#' This function tests the interaction based on a statistic proposed by Piepho (1994).
#' This function reports Piepho's test statistic, an asymptotic p-value, and a Monte Carlo p-value.
#'
#' @param x numeric matrix, \eqn{a \times b} data matrix where the number of row and column is corresponding to the number of factor levels.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#'
#'
#' @return An object of the class \code{ITtest}, which is a list inducing following components:
#' \item{pvalue.exact}{The calculated exact Monte Carlo p-value.}
#' \item{pvalue.appro}{The asymptotic p-value.}
#' \item{statistic}{The value of the test statistic.}
#' \item{Nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#' @details Piepho (1994) proposed three test statistics. The third one is
#'  based on Grubbs’ (1948) type estimator of variance for the level of the row factor.
#'  This type of estimator is used in this function. Piepho (1994) proposed an asymptotic distribution of test statistic; however, a Monte Carlo method is used to calculate the p-value.
#'  The Piepho test is not applicable when the row number of the data matrix is less than three. Note that Piepho’s test is powerful for detecting interactions when the Grubbs’ type estimators of variances are heterogeneous across the levels of one factor.
#'
#' @references Piepho, H. P. (1994). On Tests for Interaction in a Nonreplicated Two-Way Layout. Australian
#' Journal of Statistics 36:363-369.
#'
#' Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#' Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#' 86(3): 469-487.
#'
#' Grubbs, F.E. (1948). On Estimating Precision of Measuring Instruments and Product Variability. Journal of the American Statistical Association 43(242): 243-264.
#'
#' @examples
#' data(MVGH)
#' Piepho.test(MVGH, nsim = 1000)
#' 
#' @export
Piepho.test <- function(x, nsim = 10000) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    DNAME <- deparse1(substitute(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    statistics <- piephoC(x, bl, tr)
    simu <- Piephosim(nsim, bl, tr)
    pieph <- mean(statistics < simu)
    df <- bl - 1
    asypieph <- 1 - pchisq(statistics, df = df)
    out <- list(
      pvalue.exact = pieph,
      pvalue.appro = asypieph,
      nsim = nsim,
      statistic = statistics,
      data.name = DNAME,
      test = "Piepho Test"
    )
  }
  structure(out, class = "ITtest")
}
