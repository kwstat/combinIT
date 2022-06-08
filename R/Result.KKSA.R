#' @importFrom utils combn 
Result.KKSA <- function(x, nsim, alpha, simu) {
  bl <- nrow(x)
  tr <- ncol(x)
  if (bl < 4) {
    str <- paste("This test is not applicable when the row number is less than four. You may use the transpose of the data matrix if the number of column is greater than three.", '\n')
  } else {
    qKKSA <- quantile(simu, prob = alpha, names = FALSE)
    Nrow <- 2:(as.integer(bl/2))
    R <- x - matrix(rowMeans(x), bl, tr) - matrix(colMeans(x), bl, tr, byrow = TRUE) + mean(x)
    sse <- sum(R ^ 2)
    count <- 0
    Minpvalue <- 2
    f2 <- rep(0, 0)
    for(i in Nrow) {
      ind <- combn(bl, i)
      Nsplit <- ncol(ind)
      if((bl / 2) == i) Nsplit <- Nsplit / 2
      for(j in 1:Nsplit) {
        count <- count + 1
        x1 <- x[ind[, j], ]
        x2 <- x[-ind[, j], ]
        rss1 <- sum((x1 - matrix(rowMeans(x1), nrow(x1), ncol(x1)) - matrix(colMeans(x1), nrow(x1), ncol(x1), byrow = TRUE) + mean(x1))^2)
        rss2 <- sum((x2 - matrix(rowMeans(x2), nrow(x2), ncol(x2)) - matrix(colMeans(x2), nrow(x2), ncol(x2), byrow = TRUE) + mean(x2))^2)
        f2[count] <- (rss1 * (bl - i - 1))/(rss2 *(i - 1))
        if(f2[count] < 1) f2[count] <- 1 / f2[count]
        ex1 <- 1 - pf(f2[count],(tr - 1) * (i - 1), (bl - i - 1) * (tr - 1)) + pf(1/f2[count], (tr - 1) * (i - 1), (bl - i - 1) * (tr - 1))
        if(ex1 < Minpvalue) {
          Minpvalue <- ex1
          index <- ind[,j]
          RSS1 <- rss1
          RSS2 <- rss2
          df1 <- (tr - 1) * (i - 1)
          df2 <- (tr - 1) * (bl - i - 1)
          fvalue <- f2[count]
        }
      }
    }
    str1 <- paste("There may exist a significant intercation. The magnitude of interaction effects is heteroscedastic across the sub-tables of observations.", "\n")
    expre1 <- paste(index, collapse=", ")
    expre2 <- paste((1:bl)[-index], collapse =", ")
    str2 <- paste("The first sub-table consists of rows", expre1, "with RSS=", round(RSS1, 4), "on", df1, "degrees of freedoms.", "\n")
    str3 <- paste("The second sub-table consists of rows", expre2, "with RSS=", round(RSS2, 4), "on", df2, "degrees of freedoms.","\n" )
    str4 <- paste("The estimated critical value of the KKSA.test with", nsim, "Monte Carlo samples is:", round(qKKSA, 4), '\n')
    str <- paste(str1, str2, str3, str4, '\n')
  }
  str
}

