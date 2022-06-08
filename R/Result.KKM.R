<<<<<<< HEAD
=======
#' @import MASS
#' @import Matrix
>>>>>>> 73e11a32303eb884bb2a8a3ef082ec8e474b4c2c
Result.KKM <- function(x, simu, nsim, alpha, nc0) {
      qKKM <- quantile(simu, prob = 1 - alpha, names = FALSE)
      bl <- nrow(x)
      tr <- ncol(x)
      n <- bl * tr
      kp <- kpr(bl, tr)
      c0 <- C0(kp, n, nc0)
      y <- c(t(x))
      Z <- abs(kp %*% y)
      S0 <- median(Z) / c0
      PSE <- median(Z[Z <= 5 * S0])
      SZ <- Z[Z > qKKM * PSE]
      Index <- (1:nrow(kp))[Z > (qKKM * PSE)]
      if (length(Index) != 0) {
        M <- matrix(0, length(Index), 4)
        count <- 0
        for (k in Index) {
          count <- count + 1
          count2 <- 0
          for (i in 1:tr) {
            for (j in 1:bl) {
              jj <- (i - 1) * bl + j
              if (kp[k, jj] != 0) {
                count2 <- count2 + 1
                M[count, count2] <- paste0(j, i)
              }
            }
          }
        }
        C1 <- kp[Index, ]
        requireNamespace(MASS)
        requireNamespace(Matrix)
        sigma2hat <- t(y) %*% (MASS::ginv(C1) %*% C1) %*% y / Matrix::rankMatrix(C1)[1]
        str1 <- paste("There may exist a significant intercation and it might be caused by some cells.", "\n", "The absolute estimates of the significant pairwise interaction contrasts (PIC) and the corresponding involved cell means are:", "\n")
        ex1 <- paste(paste0("|mu_{", M[1, 1], "}-mu_{", M[1, 2], "}-mu_{", M[1, 3], "}+mu_{", M[1, 4], "}|="), round(SZ[1], 4), "\n")
        for (i in 2:length(Index)) {
          ex1 <- paste(ex1, paste0("|mu_{", M[i, 1], "}-mu_{", M[i, 2], "}-mu_{", M[i, 3], "}+mu_{", M[i, 4], "}|="), round(SZ[i], 4), "\n")
        }
        str2 <- ex1
        str3 <- paste("The variance estimate under the non-additivity assumption is", round(sigma2hat, 4), "on", Matrix::rankMatrix(C1)[1], "degrees of freedom.", "The estimated critical value of the KKM.test with", nsim, "Monte Carlo samples is:", round(qKKM, 4), "\n")
        str <- paste(str1, str2, str3, "\n")
      } else {
        str <- paste("The KKM.test could not detect any significant interaction.", "The estimated critical value of the KKM.test with", nsim, "Monte Carlo samples is:", round(qKKM, 4), "\n")
      }
      return(str)
}
