
Result.Piepho <- function(x, nsim, alpha, simu) {
  bl <- nrow(x)
  tr <- ncol(x)
  if (bl < 3) {
    str <- paste("This test is not applicable when the row number is less than three. You may use the transpose of the data matrix if the number of column is greater than two.", "\n")
  } else {
    qPiepho <- quantile(simu, prob = 1 - alpha, names = FALSE)
    R <- x - matrix(rowMeans(x), bl, tr) - matrix(colMeans(x), bl, tr, byrow = TRUE) + mean(x)
    W <- rowSums(R^2)
    sigmahat <- (bl * (bl - 1) * W - sum(W)) / ((bl - 1) * (bl - 2) * (tr - 1))
    str1 <- paste("There may exist a significant intercation.", "\n")
    str2 <- paste("The Grubbs' estimtors of the row variances are heterogeneous.", "\n")
    Grubbs <- paste(round(sigmahat, 4), collapse = ", ")
    str3 <- paste("The Grubbs' variance estimators are:", Grubbs, "\n")
    str4 <- paste("The estimated critical value of the Piepho.test with", nsim, "Monte Carlo samples is:", round(qPiepho, 4), "\n")
    str <- paste(str1, str2, str3, str4, "\n")
  }
  str
}
