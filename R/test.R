#' Boik (1993) locally best invariant (LBI) test
#'
#' This function calculates the LBI test statistic for testing the null hypothesis \eqn{H_0:} there is no interaction.
#' It returns an exact Monte Carlo p-value (when \eqn{p>2}) and an asymptotic chi-squared p-value.
#' 
#' @param \eqn{x} a numeric matrix, \eqn{b \times a} data matrix where the number of row and column are corresponding to the number of block and treatment levels, respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for calculating an exact Monte Carlo p-value. The default value is 10000.
#' 
#' @return  A list of consisting of:
#' @return exact.pvalue, an exact Monte Carlo p-value when \eqn{p>2}. For \eqn{p=2} an exact p-value is calculated.
#' @return asy.pvalue, a chi-squared asymptotic p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#' 
#' @details The LBI test statistic is \eqn{T_B93=(tr(R'R))^2/(p tr((R'R)^2))} where \eqn{p=min{a-1,b-1}} and \eqn{R} is the residual
#'   matrix of the input data matrix, \eqn{x}, under the null hypothesis \eqn{H_0:} there is no interaction. This test rejects the null hypothesis of no interaction when \eqn{T_B93} is small.
#'   Boik (1993) provided the exact distribution of \eqn{T_B93} when \eqn{p=2} under \eqn{H_0}. In addition, he provided an asymptotic approximation of \eqn{T_B93} under \eqn{H_0} when \eqn{q} tends to infinity where \eqn{q=max{a-1,b-1}}.
#'   Note that the LBI test is powerful when the \eqn{a \times b} matrix of interaction terms has small rank and one singular value dominates the remaining singular values or
#'   in practice, if the largest eigenvalue of \eqn{RR'} is expected to dominate the remaining eigenvalues.
#'   
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' 
#' @references Boik, R.J. (1993). Testing additivity in two-way classifications
#'  with no replications: the locally best invariant test. Journal of Applied
#'  Statistics 20(1): 41-55.
#'  
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#' @examples \dontrun{this is an example}
#' data(MVGH)
#' Boik.test(MVGH, nsim=10000)
#' @export
Boik.test <- function(x, nsim = 10000, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)
    statistics <- Bfc(x, bl, tr, p)
    simu <- Bfsim(nsim, bl, tr, p)
    boik.p <- mean(statistics > simu)
    Tb <- (1 / statistics - 1)
    T <- p * q * Tb / 2
    df <- (p + 2) * (p - 1) / 2
    if (p == 2) {
      asyboik.p <- 1 - pbeta(Tb, 1, (q - 1) / 2)
    } else {
      asyboik.p <- 1 - pchisq(T, df)
    }
    out <- list(
      exact.pvalue = boik.p, asy.pvalue = asyboik.p,
      nsim = nsim,
      statistic = statistics
    )
  }
  return(out)
}

#' Malik (2016) et al. test for interaction
#'
#' The Malik's (2016) et al. test statistics is calculated and the corresponding exact p-value is calculated by a Monte Carlo simulation.
#' 
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels, respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' 
#' @return A list of consisting of:
#' @return pvalue, the calculated exact Monte Carlo p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#' 
#' @details 
#'  Malik (2016) et al. proposed to partition
#'  the residuals into three clusters using a suitable clustering method like “k-means clustering”.
#'  The hypothesis of no interaction can be interpreted as the effect of the three
#'  clusters are equal.
#'  Note that the Malik's et al. test performs well when there are some outliers in the residuals; i.e. some cells produce large negative or positive residuals due to the significant interaction.
#'  Further, the distribution of the Malik's et al. test statistic is not known under additivity and the corresponding p-value is calculated by a Monte Carlo simulation.
#' 
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' 
#' @references Malik, W.A., Mohring, J., Piepho, H.P. (2016). A
#' clustering-based test for non-additivity in an unreplicated two-way layout.
#' Communications in Statistics-Simulation and Computation 45(2):660-670.
#' 
#' Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#' Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#' 86(3): 469-487.
#' @examples \dontrun{this is an example}
#' data(IDCP)
#' Malik.test(IDCP,nsim=10000)
#' @export
Malik.test <- function(x, nsim = 10000) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    block <- gl(bl, tr)
    treatment <- gl(tr, 1, bl * tr)
    y <- c(t(x))
    statistic <- M_f(x)
    simu <- rep(0, 0)
    for (i in 1:nsim) {
      y <- rnorm(n)
      simu[i] <- M_f(matrix(y, nrow = bl))
      # cat(paste(round(i / nsim * 100), '% completed'))
      # Sys.sleep(.05)
      # if (i == nsim) cat(': Done')
      # else cat('\014')
    }
    malik <- mean(statistic < simu)

    list(pvalue = malik, nsim = nsim, statistic = statistic)
  }
}

#' Kharrati-Kopaei and Miller (2016) test for interaction
#'
#' This function calculates the test statistic for testing \eqn{H_0:} no interaction and corresponding Monte Carlo p-value
#' proposed by Kharrati-Kopaei and Miller(2016).
#' 
#' @param \eqn{x} numerix matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels, respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param nc0 a numeric value, the number of Monte Carlo samples for computing the unbiasing constant \eqn{c_0}. The default value is 10000.
#' 
#' @return A list of consisting of:
#' @return pvalue, the calculated exact Monte Carlo p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' 
#' @details 
#' Kharrati-Kopaei and Miller(2016) proposed a test statistic for testing interaction
#' based on inspecting all pairwise interaction contrasts (PIC).
#' This test depends on an unbiasing constant \eqn{c_0} that is calculated by a Monte Carlo simulation.   
#' In addition, the null distribution of the test statistic is calculated by a Monte Carlo simulation. 
#' Note that this test procedure is powerful when significant interactions are caused by some data cells. 
#' 
#' @references Kharrati-Kopaei, M., Miller, A. (2016). A method for testing interaction in
#'  unreplicated two-way tables: using all pairwise interaction contrasts. Statistical
#'  Computation and Simulation 86(6):1203-1215.
#'   
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.

#' @examples \dontrun{this is an example}
#' data(RDWW)
#' KKM.test(RDWW,nsim=10000,nc0=10000)
#' @export
KKM.test <- function(x, nsim = 1000, nc0 = 10000, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    kp <- kpr(bl, tr)
    c0 <- C0(kp, n, nc0)
    statistics <- picf(y, kp, c0)
    simu <- PICfsim(nsim, kp, c0, n)
    PIC <- mean(statistics < simu)
    list(pvalue = PIC, nsim = nsim, statistic = statistics)
  }
}


#' Piepho (1994) test for interaction
#'
#' This function tests the interaction based on a statistic proposed by Piepho (1994).
#' This function reports Piepho's test statistic, and an asymptotic and Monte Carlo p-values.
#' 
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
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
#' @examples \dontrun{this is an example}
#' data(MVGH)
#' Piepho.test(MVGH,sim=1000)
#' @export
Piepho.test <- function(x, nsim = 10000, ...) {
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
      exact.pvalue = pieph, asy.pvalue = asypieph,
      nsim = nsim,
      statistic = statistics
    )
    return(out)
  }
}


#' Kharrati-Kopaei and Sadooghi-Alvandi (2007) test for interaction
#'
#' This function calculates Kharrati-Kopaei and Sadooghi-Alvandi's test statistic and corresponding p-value for testing interaction.
#' 
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
#'   , respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param dist a character, if dist="sim", a Monte Carlo simulation is used for calculating exact p-value. If dist="adj", the Bonferroni-adjusted p-value is calculated. The default is "sim".
#' 
#' @details  Suppose that \eqn{b>=a} and \eqn{b>=4}. Consider the \eqn{l}-th division of the data table into two sub-tables,
#'  obtained by putting \eqn{b_1} (\eqn{2≤b_1≤b-2}) rows in the first sub-table and the remaining \eqn{b_2} rows in the second sub-table (\eqn{b_1+b_2=a}).
#'  Let RSS1 and RSS2 denote the residual sum of squares for these two sub-tables, respectively. For a particular division \eqn{l}, let \eqn{F_l=max⁡(F_l,1/F_l }
#'  where \eqn{F_l=(b_2-1)RSS1/((b_1-1)RSS2)} and let \eqn{P_l} denote the corresponding p-value.
#'  Kharrati-Kopaei and Sadooghi-Alvandi (2007) proposed their test statistic as the minimum value of \eqn{P_l} over \eqn{l=1,…,2^(b-1)-b-1} all possible divisions of the table.
#'  Note that if the rows number, \eqn{b}, of data matrix is less than the columns number, \eqn{a}, the data matrix is transposed. In addition, this method of testing requires that the data matrix has more than three
#'  rows or columns. This test procedure is powerful for detecting interaction when the magnitude of interaction effects is heteroscedastic across the sub-tables of observations.
#'  
#' @return A list of consisting of:
#' @return pvalue, an exact Monte Carlo p-value.
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return statistic, the value of test statistic.
#' 
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' 
#' @references Kharrati-Kopaei, M., Sadooghi-Alvandi, S.M. (2007). A New Method for
#'  Testing Interaction in Unreplicated Two-Way Analysis of Variance. Communications
#'  in Statistics-Theory and Methods 36:2787–2803.
#'   
#'  Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#'   
#' @examples \dontrun{this is an example}
#' data(IDCP)
#' KKSA.test(IDCP,nsim=10000,dist = "sim")
#' @export
KKSA.test <- function(x, nsim = 10000, distr = "sim", ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
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
      stop("KKSA needs at least 4 levels of blocking factor")
    } else {
      cck <- 2^(bl - 1) - 1 - bl
      statistics <- kk_f(x)
      if (distr != "sim" && distr != "adj") distr <- "sim"

      if (distr == "sim") {
        simu <- rep(0, 0)
        for (i in 1:nsim) {
          simu[i] <- kk_f(matrix(rnorm(n), nrow = bl))
          cat(paste(round(i / nsim * 100), "% completed"))
          # Sys.sleep(.1)
          if (i == nsim) {
            cat(": Done", "\n")
          } else {
            cat("\014", "\n")
          }
        }
        KKSA.p <- mean(statistics > simu)
      } else if (distr == "adj") {
        KKSA.p <- statistics * cck
        KKSA.p <- min(1, KKSA.p)
      }
      out <- list(
        pvalue = KKSA.p,
        nsim = nsim, distr = distr,
        statistic = statistics
      )
      return(out)
    }
  }
}

#' Franck (2013) et al. test for interaction
#'
#' This function calculates Franck (2013) et al. test statistic,ACMIF, and corresponding p-value.
#' 
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels
#'   , respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param dist a character, if dist="sim", a Monte Carlo simulation is used for calculating exact p-value,
#'  and if dist="adj", the Bonferroni-adjusted p-value is calculated. The default is "sim".
#'  
#' @details Franck et al. (2013) derived a test statistic based on the “hidden additivity” structure.
#'  They defined this structure as “the levels of one factor belong in two or more groups such that within each group the effects of the two factors are additive but the groups may interact with the ungrouped factor”.
#'  To detect hidden additivity, Franck et al. (2013) divided the table of data into two sub-tables and an interaction F-test was developed.
#'  Then, they performed a search over all possible configuRDWWns of data and used the maximum of the interaction F-RDWWs as a test statistic. The hypothesis of no interaction is rejected when the maximum interaction F-RDWW is large.
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
#' @examples \dontrun{this is an example}
#' data(CNV)
#' Franck.test(CNV,nsim=1000,dist = "sim")
#' @export
Franck.test <- function(x, nsim = 1000, dist = "sim", ...) {
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
    if (bl < 3) {
      stop("hiddenf needs at least 3 levels of blocking factor")
    } else {
      cch <- 2^(bl - 1) - 1
      statistics <- hh_f(x)
      if (dist != "sim" & dist != "adj")
        stop("\"dist\" parameter should be equal to \"sim\" or \"adj\".")

      if (dist == "sim") {
        simu <- rep(0, 0)
        for (i in 1:nsim) {
          simu[i] <- hh_f(matrix(rnorm(n), nrow = bl))
          cat(paste(round(i / nsim * 100), "% completed", "\n"))
          # Sys.sleep(.1)
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
        nsim = nsim, dist = dist,
        statistic = statistics
      )
      return(out)
    }
  }
}

#' Combined p-value interaction test
#'
#' This function reports the p-values of the tests for non-additivity developed by Boik (1993), Piepho (1994),
#' Kharrati-Kopaei and Sadooghi-Alvandi (2007), Franck et al. (2013), Malik et al. (2016)
#' and Kharrati-Kopaei and Miller (2016). In addition, it combines the p-values of these six methods into a single p-value as a test statistic for testing interaction.
#' There are four combination methods:
#' Bonferroni, Sidak, Jacobi expansion, and Gaussian Copula. The results of these four combinations are also reported. 
#' 
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels, respectively.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param nc0 a numeric value, the number of Monte Carlo samples for computing the unbiasing constant \eqn{c_0}. The default value is 10000.
#' 
#' @details If rows number ,\eqn{b} of data matrix is less than it's columns number, \eqn{a}, 
#'  the data matrix is transposed. In addition, this test procedure requires that the data matrix has more than two
#'  rows or columns. This function Needs "mvtnorm" package.

#'  
#' @return A list of consisting of:
#' @return nsim, the number of Monte Carlo samples that are used to estimate p-value.
#' @return piepho.pvalue, the p-value of Piepho's (1994) test.
#' @return Boik.pvalue, the p-value of Boik's (1993) test.
#' @return Malik.pvalue, the p-value of Malik's (2016) et al. test.
#' @return KKM.pvalue, the p-value of Kharrati-Kopaei and Miller's (2016) test.
#' @return KKSA.pvalue, the p-value of Kharrati-Kopaei and Sadooghi-Alvandi's (2007) test.
#' @return Franck.pvalue, the p-value of Franck's (2013) et al. test.
#' @return Bonferroni, the combined p-value by using the Bonferroni method.
#' @return Sidak, the combined p-value by using the Sidak method.
#' @return jacobi, the combined p-value by using the Jacobi method.
#' @return GC, the combined p-value by using the Guassian copula.
#' 
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' 
#' @references Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#' @examples \dontrun{this is an example}
#' data(RDWW)
#' CPI.test(RDWW,nsim=500,nc0=10000)
#' @export
CPI.test <- function(x, nsim = 500, nc0 = 10000, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    if (bl < tr) warning("transpose the input matrix")
    x <- t(x)
    te <- bl
    bl <- tr
    tr <- te
  }
  if (bl < 3) {
    stop("hiddenf needs at least 3 levels of blocking factor")
  } else {
    n <- tr * bl
    block <- gl(bl, tr)
    treatment <- gl(tr, 1, bl * tr)
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)
    cck <- 2^(bl - 1) - 1 - bl
    cch <- 2^(bl - 1) - 1
    kp <- kpr(bl, tr)
    c0 <- mean(replicate(nc0, {
      median(abs(kp %*% rnorm(n)))
    }))

    sta <- bmp.f(x, y, block, treatment, bl, tr, p)

    Bstat <- sta$Boik
    Mstat <- sta$Tc
    pistat <- sta$piepho
    pstat <- pic.f(y, kp, c0)
    if (bl == 3) {
      Hstat <- hh_f(x, bl)
    } else {
      Ksimu <- rep(0, 0)
      kh <- kh_f(x)
      Kstat <- kh$fmin
      Hstat <- kh$fmax
    }

    Bsimu <- Msimu <- psimu <- pisimu <- Hsimu <- rep(0, 0)
    for (i in 1:nsim) {
      y <- rnorm(n)
      x <- matrix(y, nrow = bl, byrow = TRUE)
      sta <- bmp.f(x, y, block, treatment, bl, tr, p)
      Bsimu[i] <- sta$Boik
      Msimu[i] <- sta$Tc
      pisimu[i] <- sta$piepho
      psimu[i] <- pic.f(y, kp, c0)
      if (bl == 3) {
        Hsimu[i] <- hh_f(x, bl)
      } else {
        kh <- kh_f(x)
        Ksimu[i] <- kh$fmin
        Hsimu[i] <- kh$fmax
      }
      cat(paste(round(i / nsim * 100), "% completed"), "\n")
      Sys.sleep(.1)
      if (i == nsim) {
        cat(": Done", "\n")
      } else {
        cat("\014", "\n")
      }
    }
    Boik.pvalue <- mean(Bstat > Bsimu)
    piepho.pvalue <- mean(pistat < pisimu)
    PIC.pvalue <- mean(pstat < psimu)
    Malik.pvalue <- mean(Mstat < Msimu)
    hiddenf.pvalue <- mean(Hstat < Hsimu)
    if (bl == 3) {
      KKSA.pvalue <- NULL
    } else {
      KKSA.pvalue <- mean(Kstat > Ksimu)
    }
    pvalues <- c(Boik.pvalue, piepho.pvalue, hiddenf.pvalue, Malik.pvalue, PIC.pvalue, KKSA.pvalue)
    cp <- comb(pvalues)
    Bonferroni <- cp$Bon
    GC <- cp$GC
    Sidak <- cp$Sidak
    jacobi <- cp$jacobi
    list(
      nsim = nsim, piepho.pvalue = piepho.pvalue, Boik.pvalue = Boik.pvalue,
      Malik.pvalue = Malik.pvalue, KKM.pvalue = PIC.pvalue,
      KKSA.pvalue = KKSA.pvalue, Franck.pvalue = hiddenf.pvalue,
      Bonferroni = Bonferroni, Sidak = Sidak, jacobi = jacobi, GC = GC
    )
  }
}


#' Interaction plot
#'
#' @param \eqn{x} numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels, respectively.
#' @return  An interaction plot for input
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' @examples \dontrun{this is an example}
#' data(CNV)
#' interactionplot(CNV)
#' @export
interactionplot <- function(x, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    par(mfcol = c(1, 2))
    t <- ncol(x)
    b <- nrow(x)
    matplot(t(x), type = "b", xaxt = "n", ylab = "Observed values", xlab = "Factor1(column)", lty = 1:b, ...)
    axis(1, at = 1:t, labels = 1:t, cex.axis = 1)
    legend("topright", rep(paste0("row", 1:b)), lty = 1:b, bty = "n", cex = 0.7)
    matplot(x, type = "b", xaxt = "n", ylab = "", xlab = "Factor2(row)", lty = 1:t, ...)
    legend("topright", rep(paste0("col", 1:t)), lty = 1:t, bty = "n", cex = 0.7)
    axis(1, at = 1:b, labels = 1:b, cex.axis = 1)
  }
}