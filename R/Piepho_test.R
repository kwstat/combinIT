#' Piepho (1994) Test for Interaction
#'
#' This function tests the interaction based on a statistic proposed by Piepho (1994).
#' This function reports Piepho's test statistic, and an asymptotic and Monte Carlo p-values.
#'
#' @param x numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
#'   , respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#'
#' @return A list of consisting of:
#' @return exact.pvalue, an exact Monte Carlo p-value.
#' @return asy.pvalue, an asymptotic p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#'
#' @details Piepho (1994) proposed three test statistics.The third one is
#'  based on Grubbs’ (1948) type estimator of variance for each level of block effect.
#'  This type of estimator is used in this function. Piepho (1994) proposed an asymptotic distribution of test statistic; however, we use a Monte Carlo method to calculate the p-value.
#'  Note that Piepho’s test is powerful for detecting interactions when the Grubbs’ type estimators of variances are heterogeneous across the levels of one factor.
#'
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
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
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    statistics <- piephoC(x, bl, tr)
    simu <- Piephosim(nsim, bl, tr)
    pieph <- mean(statistics < simu)
    df <- bl - 1
    asypieph <- 1 - pchisq(statistics, df = df)
    out <- list(
      exact.pvalue = pieph,
      asy.pvalue = asypieph,
      nsim = nsim,
      statistic = statistics
    )
    return(out)
  }
}
