#' Combined P-value Interaction Test
#'
#' This function reports the p-values of the tests for non-additivity developed by Boik (1993), Piepho (1994),
#' Kharrati-Kopaei and Sadooghi-Alvandi (2007), Franck et al. (2013), Malik et al. (2016)
#' and Kharrati-Kopaei and Miller (2016). In addition, it combines the p-values of these six tests (and some other available p-values) into a single p-value as a test statistic for testing interaction.
#' There are four combination methods:
#' Bonferroni, Sidak, Jacobi expansion, and Gaussian Copula. The results of these four combined tests are also reported. If there is a significant interaction, the type of interaction is also provided.
#'
#' @param x numeric matrix, \eqn{a \times b} data matrix where the number of row and column is corresponding to the number of factor levels.
#' @param nsim a numeric value, the number of Monte Carlo samples for computing an exact Monte Carlo p-value. The default value is 10000.
#' @param nc0 a numeric value, the number of Monte Carlo samples for computing the unbiasing constant \eqn{c_0} in \code{KKM.test}. The default value is 10000.
#' @param opvalue a numeric vector, other p-values (in addition to the six considered p-values) that are going to be combined.
#' @param alpha a numeric value, the level of the test. The default value is 0.05.
#' @param report logical: if \code{TRUE} the result of the test is reported at the \code{alpha} level.
#' @param Elapsed.time logical: if \code{TRUE} the progress will be printed in the console.
#'
#' @details The data matrix is divided based on the row of the data matrix for \code{KKSA.test} and \code{Franck.test}. Note that \code{KKSA.test} is not applicable when \eqn{a} is less than four. \code{Franck.test} and \code{Piepho.test} are not applicable when \eqn{a} is less than three. This function needs \code{mvtnorm} package.
#'
#' @return An object of the class \code{combtest}, which is a list inducing following components::
#' \item{nsim}{The number of Monte Carlo samples that are used to estimate p-value.}
#' \item{Piepho.pvalue}{The p-value of Piepho's (1994) test.}
#' \item{Piepho.Stat}{The value of Piepho's (1994) test statistic.}
#' \item{Boik.pvalue}{The p-value of Boik's (1993) test.}
#' \item{Boik.Stat}{The value of Boik's (1993) test statistic.}
#' \item{Malik.pvalue}{The p-value of Malik's (2016) et al. test.}
#' \item{Malik.Stat}{The value of Malik's (2016) et al. test statistic.}
#' \item{KKM.pvalue}{The p-value of Kharrati-Kopaei and Miller's (2016) test.}
#' \item{KKM.Stat}{The value of Kharrati-Kopaei and Miller's (2016) test statistic.}
#' \item{KKSA.pvalue}{The p-value of Kharrati-Kopaei and Sadooghi-Alvandi's (2007) test.}
#' \item{KKSA.Stat}{The value of Kharrati-Kopaei and Sadooghi-Alvandi's (2007) test statistic.}
#' \item{Franck.pvalue}{The p-value of Franck's (2013) et al. test.}
#' \item{Franck.Stat}{The value of Franck's (2013) et al. test statistic.}
#' \item{Bonferroni}{The combined p-value by using the Bonferroni method.}
#' \item{Sidak}{The combined p-value by using the Sidak method.}
#' \item{Jacobi}{The combined p-value by using the Jacobi method.}
#' \item{GC}{The combined p-value by using the Gaussian copula.}
#' \item{data.name}{The name of the input dataset.}
#' \item{test}{The name of the test.}
#' \item{Level}{The level of test.}
#' \item{Result}{The result of the combined test at the alpha level with some descriptions on the type of significant interaction.}
#'
#' @references Shenavari, Z., Kharrati-Kopaei, M. (2018). A Method for Testing Additivity in
#'  Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests. International Statistical Review
#'  86(3): 469-487.
#'
#' @examples
#' data(RDWW)
#' CPI.test(RDWW, nsim = 1000, Elapsed.time = FALSE)
#'
#' @importFrom stats pchisq pf qnorm var
#'
#' @export
CPI.test <- function(x, nsim = 10000, nc0 = 10000, opvalue = NULL, alpha = 0.05, report = TRUE, Elapsed.time = TRUE) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    DNAME <- deparse1(substitute(x))
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    if (bl == 3) {
      warning("KKSA.test needs at least 4 levels for the row factor. For combining pvalues, the pvalue of the KKSA test is not considered.")
    }
    if (bl <= 2) {
      warning("Franck.test and Piepho.test need at least 3 levels for the row factor. KKSA.test also needs at least 4 levels for the row factor. For combining pvalues, the pvalues of the Franck, Piepho, and KKSA tests are not considered.")
    }
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
    sta <- bmp_f(x)
    Bstat <- sta$Boik
    Mstat <- sta$Tc
    pistat <- sta$piepho
    pstat <- picf(y, kp, c0)
    if (bl == 3) {
      Hstat <- hh_f(x)
    }
    if (bl > 3) {
      Ksimu <- rep(0, 0)
      kh <- kh_f(x)
      Kstat <- kh$fmin
      Hstat <- kh$fmax
    }
    Bsimu <- Msimu <- psimu <- pisimu <- Hsimu <- rep(0, 0)
    if (Elapsed.time) {
      pb <- completed(nsim)
      for (i in 1:nsim) {
        yy <- rnorm(n)
        xx <- matrix(yy, nrow = bl, byrow = TRUE)
        sta <- bmp_f(xx)
        Bsimu[i] <- sta$Boik
        Msimu[i] <- sta$Tc
        pisimu[i] <- sta$piepho
        psimu[i] <- picf(yy, kp, c0)
        if (bl == 3) {
          Hsimu[i] <- hh_f(xx)
        }
        if (bl > 3) {
          kh <- kh_f(xx)
          Ksimu[i] <- kh$fmin
          Hsimu[i] <- kh$fmax
        }
        if (i == pb$pr[pb$j]) pb <- nextc(pb, i)
      }
    } else {
      for (i in 1:nsim) {
        yy <- rnorm(n)
        xx <- matrix(yy, nrow = bl, byrow = TRUE)
        sta <- bmp_f(xx)
        Bsimu[i] <- sta$Boik
        Msimu[i] <- sta$Tc
        pisimu[i] <- sta$piepho
        psimu[i] <- picf(yy, kp, c0)
        if (bl == 3) {
          Hsimu[i] <- hh_f(xx)
        }
        if (bl > 3) {
          kh <- kh_f(xx)
          Ksimu[i] <- kh$fmin
          Hsimu[i] <- kh$fmax
        }
      }
    }
    Tb <- (1 / Bstat - 1)
    if (p == 1) {
      Boik.pvalue <- 1
    }
    if (p == 2) {
      Boik.pvalue <- 1 - pbeta(Tb, 1, (q - 1) / 2)
    }
    if (p > 2) {
      Boik.pvalue <- mean(Bstat >= Bsimu)
    }
    PIC.pvalue <- mean(pstat < psimu)
    Malik.pvalue <- mean(Mstat < Msimu)
    qMalik <- quantile(Msimu, prob = 1 - alpha, names = FALSE)
    qBoik <- quantile(Bsimu, prob = alpha, names = FALSE)
    if (bl <= 3) {
      KKSA.pvalue <- NA
    } else {
      KKSA.pvalue <- mean(Kstat > Ksimu)
      qKKSA <- quantile(Ksimu, prob = alpha, names = FALSE)
    }
    if (bl <= 2) {
      hiddenf.pvalue <- NA
      piepho.pvalue <- NA
    } else {
      hiddenf.pvalue <- mean(Hstat < Hsimu)
      piepho.pvalue <- mean(pistat < pisimu)
      qPiepho <- quantile(pisimu, prob = 1 - alpha, names = FALSE)
      qFranck <- quantile(Hsimu, prob = 1 - alpha, names = FALSE)
    }
    pvalues <- c(Boik.pvalue, piepho.pvalue, hiddenf.pvalue, Malik.pvalue, PIC.pvalue, KKSA.pvalue, opvalue)
    if (is.null(opvalue)) {
      names(pvalues) <- c("Boik.test", "Piepho.test", "Franck.test", "Malik.test", "KKM.test", "KKSA.test")
    } else {
      names(pvalues) <- c("Boik.test", "Piepho.test", "Franck.test", "Malik.test", "KKM.test", "KKSA.test", paste0("Added test", 1:length(opvalue)))
    }
    if (bl <= 3) {
      pvalues <- pvalues[!is.na(pvalues)]
    } else {
      pvalues <- pvalues
    }
    spvalues <- sort(pvalues)
    cp <- comb(pvalues)
    Bonferroni <- cp$Bon
    GC <- cp$GC
    Sidak <- cp$Sidak
    jacobi <- cp$jacobi
    if (report) {
      if (cp$Bon >= alpha & cp$GC >= alpha & cp$Sidak >= alpha & cp$jacobi >= alpha) {
        str <- paste("No significant interaction was detected at the", paste0(100 * alpha, "%"), "level.", "\n")
      }
      if ((cp$Bon < alpha | cp$Sidak < alpha | cp$jacobi < alpha) & bl >= 4) {
        if (min(pvalues) == Boik.pvalue) str <- Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu)
        if (min(pvalues) == piepho.pvalue) str <- Result.Piepho(x, nsim = nsim, alpha = alpha, simu = pisimu)
        if (min(pvalues) == hiddenf.pvalue) str <- Result.Franck(x, nsim = nsim, alpha = alpha, simu = Hsimu)$string
        if (min(pvalues) == Malik.pvalue) str <- Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha)
        if (min(pvalues) == PIC.pvalue) str <- Result.KKM(x, nsim = nsim, simu = psimu, alpha = alpha, nc0 = nc0)
        if (min(pvalues) == KKSA.pvalue) str <- Result.KKSA(x, nsim = nsim, alpha = alpha, simu = Ksimu)$string
        if (any(min(pvalues) == opvalue)) str <- paste("Significant interactions may be due to the", paste0(names(spvalues[1])), "that its p-value is recently added.")
      }
      if ((cp$Bon < alpha | cp$Sidak < alpha | cp$jacobi < alpha) & bl == 3) {
        if (min(pvalues) == Boik.pvalue) str <- Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu)
        if (min(pvalues) == piepho.pvalue) str <- Result.Piepho(x, nsim = nsim, alpha = alpha, simu = pisimu)
        if (min(pvalues) == hiddenf.pvalue) str <- Result.Franck(x, nsim = nsim, alpha = alpha, simu = Hsimu)$string
        if (min(pvalues) == Malik.pvalue) str <- Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha)
        if (min(pvalues) == PIC.pvalue) str <- Result.KKM(x, simu = psimu, alpha = alpha, nc0 = nc0)
        if (any(min(pvalues) == opvalue)) str <- paste("Significant interactions may be due to the", paste0(names(spvalues[1])), "that its p-value is recently added.")
      }
      if ((cp$Bon < alpha | cp$Sidak < alpha | cp$jacobi < alpha) & bl < 3) {
        if (min(pvalues) == Boik.pvalue) str <- Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu)
        if (min(pvalues) == Malik.pvalue) str <- Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha)
        if (min(pvalues) == PIC.pvalue) str <- Result.KKM(x, simu = psimu, alpha = alpha, nc0 = nc0)
        if (any(min(pvalues) == opvalue)) str <- paste("Significant interactions may be due to the", paste0(names(spvalues[1])), "that its p-value is recently added.")
      }
      if (cp$Bon < alpha) {
        for (i in 2:length(spvalues)) {
          if (spvalues[i] < alpha / (length(spvalues)-i+1)) {
            str <- paste(str, "\n", "In addition to the", paste0(names(spvalues[i-1]), ", the"), names(spvalues[i]), "is significant", "by using the Holm-Bonferroni method:")
            if (bl >= 4) {
              if (spvalues[i] == Boik.pvalue) str <- paste(str, Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu))
              if (spvalues[i] == piepho.pvalue) str <- paste(str, Result.Piepho(x, nsim = nsim, alpha = alpha, simu = pisimu))
              if (spvalues[i] == hiddenf.pvalue) str <- paste(str, Result.Franck(x, nsim = nsim, alpha = alpha, simu = Hsimu)$string)
              if (spvalues[i] == Malik.pvalue) str <- paste(str, Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha))
              if (spvalues[i] == PIC.pvalue) str <- paste(str, Result.KKM(x, nsim = nsim, simu = psimu, alpha = alpha, nc0 = nc0))
              if (spvalues[i] == KKSA.pvalue) str <- paste(str, Result.KKSA(x, nsim = nsim, alpha = alpha, simu = Ksimu)$string)
              if (any(spvalues[i] == opvalue)) str <- paste(str, paste("Significant interactions may be due to the", paste0(names(spvalues[i])), "that its p-value is recently added."))
            }
            if (bl == 3) {
              if (spvalues[i] == Boik.pvalue) str <- paste(str, Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu))
              if (spvalues[i] == piepho.pvalue) str <- paste(str, Result.Piepho(x, nsim = nsim, alpha = alpha, simu = pisimu))
              if (spvalues[i] == hiddenf.pvalue) str <- paste(str, Result.Franck(x, nsim = nsim, alpha = alpha, simu = Hsimu)$string)
              if (spvalues[i] == Malik.pvalue) str <- paste(str, Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha))
              if (spvalues[i] == PIC.pvalue) str <- paste(str, Result.KKM(x, simu = psimu, alpha = alpha, nc0 = nc0))
              if (any(spvalues[i] == opvalue)) str <- paste(str, paste("Significant interactions may be due to the", paste0(names(spvalues[i])), "that its p-value is recently added."))
            }
            if (bl < 3) {
              if (spvalues[i] == Boik.pvalue) str <- paste(str, Result.Boik(x, nsim = nsim, alpha = alpha, simu = Bsimu))
              if (spvalues[i] == Malik.pvalue) str <- paste(str, Result.Malik(x, simu = Msimu, nsim = nsim, alpha = alpha))
              if (spvalues[i] == PIC.pvalue) str <- paste(str, Result.KKM(x, simu = psimu, alpha = alpha, nc0 = nc0))
              if (any(spvalues[i] == opvalue)) str <- paste(str, paste("Significant interactions may be due to the", paste0(names(spvalues[i])), "that its p-value is recently added."))
            }
          }
        }
      }
    } else {
      str <- paste("A report has not been wanted! To have a report, change argument 'report' to TRUE.")
    }
    if (bl < 4) {
      KKSA.pvalue <- NA
      Kstat <- NA
    }
    if (bl < 3) {
      hiddenf.pvalue <- NA
      Hstat <- NA
      piepho.pvalue <- NA
      pistat <- NA
    }
    out <- list(
      nsim = nsim,
      Piepho.pvalue = piepho.pvalue,
      Piepho.Stat = pistat,
      Boik.pvalue = Boik.pvalue,
      Boik.Stat = Bstat,
      Malik.pvalue = Malik.pvalue,
      Malik.Stat = Mstat,
      KKM.pvalue = PIC.pvalue,
      KKM.Stat = pstat,
      KKSA.pvalue = KKSA.pvalue,
      KKSA.Stat = Kstat,
      Franck.pvalue = hiddenf.pvalue,
      Franck.Stat = Hstat,
      Bonferroni = Bonferroni,
      Sidak = Sidak,
      Jacobi = jacobi,
      GC = GC,
      data.name = DNAME,
      test = "Combined p-value interaction Test",
      Level = alpha,
      Result = str
    )
    structure(out, class = "combtest")
  }
}
