
Result.Malik <- function(x, simu, alpha, nsim) {
  bl <- nrow(x)
  tr <- ncol(x)
  qMalik <- quantile(simu, prob = 1 - alpha, names = FALSE)
  R <- x - matrix(rowMeans(x), bl, tr) - matrix(colMeans(x), bl, tr, byrow = TRUE) + mean(x)
  R <- round(R, 5)
  Min <- min(R)
  Max <- max(R)
  Index.Min <- which(R == Min)
  Index.Max <- which(R == Max)
  pmax <- Index.Max %/% bl
  rmax <- Index.Max %% bl
  pmin <- Index.Min %/% bl
  rmin <- Index.Min %% bl
  cellmax <- matrix(0, nrow = length(pmax), ncol = 2)
  cellmin <- matrix(0, nrow = length(pmin), ncol = 2)
  for (i in 1:length(pmax)) {
    if (rmax[i] == 0) {
      cellmax[i, ] <- c(bl, pmax[i])
    } else {
      cellmax[i, ] <- c(rmax[i], pmax[i] + 1)
    }
  }
  for (i in 1:length(pmin)) {
    if (rmin[i] == 0) {
      cellmin[i, ] <- c(bl, pmin[i])
    } else {
      cellmin[i, ] <- c(rmin[i], pmin[i] + 1)
    }
  }
  str1 <- paste("There may exist a significant intercation.", "The significant interaction might due to the some outliers in residuals; some cells produce large negative or positive residuals.", "\n")
  str2 <- paste("The cell with row=", cellmin[1, 1], "and column=", cellmin[1, 2], "produces a large negative residual.", "\n")
  if (length(pmin) > 2) {
    for (i in 2:length(pmin)) {
      str2 <- paste(str2, "The cell with row=", cellmin[i, 1], "and column=", cellmin[i, 2], "produces a large negative residual.", "\n")
    }
  }
  str3 <- paste("The cell with row=", cellmax[1, 1], "and column=", cellmax[1, 2], "produces a large positive residual.", "\n")
  if (length(pmax) > 2) {
    for (i in 2:length(pmax)) {
      str3 <- paste(str3, "The cell with row=", cellmax[i, 1], "and column=", cellmax[i, 2], "produces a large positive residual.", "\n")
    }
  }
  str4 <- paste("The estimated critical value of the Malik.test with", nsim, "Monte Carlo samples is:", round(qMalik, 4), "\n")
  str <- paste(str1, str2, str3, str4, "\n")
  str
}
