#' @importFrom utils combn
Result.Franck <- function(x, nsim, alpha, simu) {
  tr <- ncol(x)
  bl <- nrow(x)
  if (bl < 3) {
    str <- paste("This test is not applicable when the row number is less than three. You may use the transpose of the data matrix if the number of column is greater than two.", "\n")
    sb1 <- 1:bl
  } else {
    qFranck <- quantile(simu, prob = 1 - alpha, names = FALSE)
    cc <- 2^(bl - 1) - 1
    Nrow <- 2:(as.integer(bl / 2))
    sse <- sum((x - matrix(rowMeans(x), bl, tr) - matrix(colMeans(x), bl, tr, byrow = TRUE) + mean(x))^2)
    fvalues <- rep(0, cc)
    count <- 0
    maxfv <- 0
    for (i in Nrow) {
      ind <- combn(bl, i)
      Nsplit <- ncol(ind)
      if (bl / 2 == i) Nsplit <- Nsplit / 2
      for (j in 1:Nsplit) {
        count <- count + 1
        x1 <- x[ind[, j], ]
        if (length(ind[, j]) == 1){
          x1 <- matrix(x1, 1, ncol(x))
        }
        x2 <- x[(1:bl)[-ind[, j]], ]
        if (length((1:bl)[-ind[, j]]) == 1) {
          x2 <- matrix(x2, 1, ncol(x))
        }
        rss1 <- sum((x1 - matrix(rowMeans(x1), nrow(x1), ncol(x1)) - matrix(colMeans(x1), nrow(x1), ncol(x1), byrow = TRUE) + mean(x1))^2)
        rss2 <- sum((x2 - matrix(rowMeans(x2), nrow(x2), ncol(x2)) - matrix(colMeans(x2), nrow(x2), ncol(x2), byrow = TRUE) + mean(x2))^2)
        sse7 <- rss1 + rss2
        fvalues[count] <- (sse - sse7) * (bl - 2) / sse7
        if (fvalues[count] > maxfv) {
          maxfv <- fvalues[count]
          sb1 <- ind[, j]
        }
      }
    }
    for (d in 1:bl) {
      count <- count + 1
      x3 <- x[-d, ]
      sse7 <- sum((x3 - matrix(rowMeans(x3), nrow(x3), ncol(x3)) - matrix(colMeans(x3), nrow(x3), ncol(x3), byrow = TRUE) + mean(x3))^2)
      fvalues[count] <- (sse - sse7) * (bl - 2) / sse7
      if (fvalues[count] > maxfv) {
        maxfv <- fvalues[count]
        sb1 <- c(1:bl)[-d]
      }
    }
    sb2 <- c(1:bl)[-sb1]
    str1 <- paste("A significant hidden structure of intercation might exist.")
    expre1 <- paste((sb1), collapse = ", ")
    expre2 <- paste((sb2), collapse = ", ")
    str2 <- paste("The first group includes rows:", expre1,".")
    str3 <- paste("The second group includes rows:", expre2,".")
    str4 <- paste("The estimated critical value of the Franck.test with", nsim, "Monte Carlo samples is", round(qFranck, 4),".")
    str <- paste(str1, str2, str3, str4)
  }
  list(string = str, index = sb1)
}
