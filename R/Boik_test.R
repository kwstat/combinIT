#' Boik's (1993) Locally Best Invariant (LBI) Test
#'
#' This function calculates the LBI test statistic for testing the null hypothesis \eqn{H_0:} there is no interaction.
#' It returns an exact p-value when \eqn{p=2}. It returns an exact Monte Carlo p-value when \eqn{p>2}. It also provides an asymptotic chi-squared p-value. Note that the p-value of the Boik.test is always 1 when \eqn{p=1}.
#'
#' @param x a numeric matrix, \eqn{b \times a} data matrix where the number of row and column are corresponding to the number of block and treatment levels, respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for calculating an exact Monte Carlo p-value. The default value is 10000.
#'
#' @return An object of the class \code{ITtest}, which is a list inducing following components::
#' \item{pvalue.exact}{An exact Monte Carlo p-value when \eqn{p>2}. For \eqn{p=2} an exact p-value is calculated.}
#' \item{pvalue.appro}{An chi-squared asymptotic p-value.}
#' \item{statistic}{The value of test statistic.}
#' \item{Nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#'
#'
#' @details The LBI test statistic is \eqn{T_B93=(tr(R'R))^2/(p tr((R'R)^2))} where \eqn{p=min{a-1,b-1}} and \eqn{R} is the residual
#'   matrix of the input data matrix, \eqn{x}, under the null hypothesis \eqn{H_0:} there is no interaction. This test rejects the null hypothesis of no interaction when \eqn{T_B93} is small.
#'   Boik (1993) provided the exact distribution of \eqn{T_B93} when \eqn{p=2} under \eqn{H_0}. In addition, he provided an asymptotic approximation of \eqn{T_B93} under \eqn{H_0} when \eqn{q} tends to infinity where \eqn{q=max{a-1,b-1}}.
#'   Note that the LBI test is powerful when the \eqn{a \times b} matrix of interaction terms has small rank and one singular value dominates the remaining singular values or
#'   in practice, if the largest eigenvalue of \eqn{RR'} is expected to dominate the remaining eigenvalues.
#'
#'
#' @references Boik, R.J. (1993). Testing additivity in two-way classifications
#'  with no replications: the locally best invariant test. Journal of Applied
#'  Statistics 20(1): 41-55.
#'
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#'
#' @examples
#' data(MVGH)
#' Boik.test(MVGH, nsim = 1000)
#' @importFrom stats median pbeta rnorm
#' @export
Boik.test <- function(x, nsim = 10000) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    DNAME <- deparse1(substitute(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)
    statistics <- Bfc(x, bl, tr, p)
    Tb <- (1 / statistics - 1)
    T0 <- p * q * Tb / 2
    df <- (p + 2) * (p - 1) / 2
    if (p == 1) {
      asyboik.p <- 1
      simu <- Bfsim(nsim, bl, tr, p)
      boik.p <- mean(statistics >= simu)
    }
    if (p > 2) {
      simu <- Bfsim(nsim, bl, tr, p)
      boik.p <- mean(statistics >= simu)
      asyboik.p <- 1 - pchisq(T0, df)
    }
    if (p == 2) {
      boik.p <- 1 - pbeta(Tb, 1, (q - 1) / 2)
      asyboik.p <- 1 - pchisq(T0, df)
    }
    out <- list(
      pvalue.exact = boik.p,
      pvalue.appro = asyboik.p,
      nsim = nsim,
      statistic = statistics,
      data.name = DNAME,
      test = "Boik Test"
    )
  }
  structure(out, class = "ITtest")
}
