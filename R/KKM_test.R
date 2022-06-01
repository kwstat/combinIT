#' Kharrati-Kopaei and Miller's (2016) Test for Interaction
#'
#' This function calculates the test statistic for testing \eqn{H_0:} There is no interaction, and corresponding Monte Carlo p-value
#' proposed by Kharrati-Kopaei and Miller(2016).
#'
#' @param x a numeric matrix, \eqn{a \times b} data matrix where the number of row and column is corresponding to the number of factor levels.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param nc0 a numeric value, the number of Monte Carlo samples for computing the unbiasing constant \eqn{c_0}. The default value is 10000.
#' @param alpha a numeric value, the level of the test. The default value is 0.05.
#' 
#'
#' @return An object of the class \code{ITtest}, which is a list inducing following components::
#' \item{pvalue.exact}{The calculated exact Monte Carlo p-value.}
#' \item{pvalue.appro}{is not available for \code{KKM.test}.}
#' \item{Nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{statistic}{The value of the test statistic.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#' \item{Level}{The level of test.}
#' \item{Result}{The result of the test at the alpha level with some descriptions on the type of significant interaction.}
#' 
#' @details
#' Kharrati-Kopaei and Miller(2016) proposed a test statistic for testing interaction
#' based on inspecting all pairwise interaction contrasts (PIC).
#' This test depends on an unbiasing constant \eqn{c_0} that is calculated by a Monte Carlo simulation.
#' In addition, the null distribution of the test statistic is calculated by a Monte Carlo simulation. This test is not applicable when both \eqn{a} and \eqn{b} are less than three.
#' Note that this test procedure is powerful when significant interactions are caused by some data cells.
#'
#' @references Kharrati-Kopaei, M., Miller, A. (2016). A method for testing interaction in
#'  unreplicated two-way tables: using all pairwise interaction contrasts. Statistical
#'  Computation and Simulation 86(6):1203-1215.
#'
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.

#' @examples
#' data(RDWW)
#' KKM.test(RDWW, nsim = 1000, nc0 = 1000)
#' 
#' @export
KKM.test <- function(x, nsim = 1000, alpha = 0.05, nc0 = 10000) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    library(Matrix)
    library(MASS)
    DNAME <- deparse1(substitute(x))
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    kp <- kpr(bl, tr)
    c0 <- C0(kp, n, nc0)
    statistics <- picf(y, kp, c0)
    simu <- PICfsim(nsim, kp, c0, n)
    PIC <- mean(statistics < simu)
    if (PIC < alpha) {
      qKKM <- quantile(simu, prob = 1 - alpha, names = FALSE)
      Z <- abs(kp %*% y)
      S0 <- median(Z) / c0
      PSE <- median(Z[Z <= 5 * S0])
      SZ <- Z[Z > qKKM * PSE]
      Index <- (1:nrow(kp))[Z > (qKKM * PSE)]
      M <- matrix(0, length(Index), 4)
      count <- 0
      for(k in Index) {
        count <- count + 1
        count2 <- 0
        for (i in 1:tr) {
          for (j in 1:bl) {
            jj <- (i-1) *tr + j
            if (kp[k, jj] != 0) {
              count2 <- count2 + 1
              M[count, count2] <- paste0(j, i)
            }
          }
        }
      }
      C1 <- kp[Index, ]
      sigma2hat <- t(y) %*% (ginv(C1) %*% C1) %*% y / rankMatrix(C1)[1]
      str1 <- paste("There may exist a significant intercation and it might be caused by some cells.", '\n', "The absolute estimates of the significant pairwise interaction contrasts (PIC) and the corresponding involved cell means are:", '\n')
      ex1 <- paste("|mu_{", M[1,1], "}-mu_{", M[1,2], "}-mu_{", M[1,3], "}+mu_{", M[1,4], "}|=", round(SZ[1], 4), '\n')
      for (i in 2:length(Index)) {
          ex1 <- paste(ex1, "|mu_{", M[i,1], "}-mu_{", M[i,2], "}-mu_{", M[i,3], "}+mu_{", M[i,4], "}|=", round(SZ[i], 4), '\n')
      }
      str2 <- ex1
      str3 <- paste("The estimate of the variance under the non-additivity assumption is", round(sigma2hat, 4), "on", rankMatrix(C1)[1], "degrees of freedom", '\n')
      str <- paste(str1, str2, str3)
    } else {
      str <- "The KKM.test could not detect any significant interaction."
    } 
  }
  structure(
    list(
      pvalue.exact = PIC,
      pvalue.appro = "NULL",
      nsim = nsim,
      statistic = statistics,
      data.name = DNAME,
      test = "KKM Test",
      Level = alpha,
      Result = str
    ),
    class = "ITtest"
  )
}