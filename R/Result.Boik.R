#' @importFrom stats qbeta
Result.Boik <- function(x, nsim, alpha, simu) {
  bl <- nrow(x)
  tr <- ncol(x)
  n <- tr * bl
  p <- min(tr - 1, bl - 1)
  q <- max(tr - 1, bl - 1)
  df <- (p + 2) * (p - 1) / 2
  if (p == 1) {
    boik.p <- 1
    qBoik <- 1
    str3 <- paste("The exact critical value of the Boik.test is:", 1, '\n')
  }
  if (p > 2) {
    qBoik <- quantile(simu, prob = alpha, names = FALSE)
    str3 <- paste("The estimated critical value of the Boik.test with", nsim, "Monte Carlo samples is:", round(qBoik, 4), '\n')    
  }
  if (p == 2) {
    qBoik <- qbeta(1 - alpha, 1, (q - 1) /2)
    qBoik <- 1/(qBoik + 1)
    str3 <- paste("The exact critical value of the Boik.test is:", round(qBoik, 4), '\n')    
  }
  R <- x - matrix(rowMeans(x), bl, tr) - matrix(colMeans(x), bl, tr, byrow = TRUE) + mean(x)
  EV <- round(eigen(R%*%t(R))$values, 4)
  str1 <- paste("There may exist a significant multiplicative form of intercation.",'\n')
  expre <- paste(as.character(EV), collapse = ", ")
  str2 <- paste("The eigen values of the RR' matrix are:", expre, '\n')
  str3 <- str3 #paste("The estimated critical value of the Boik.test with", nsim, "Monte Carlo samples is:", round(qBoik, 4), '\n')
  str <- paste(str1, str2, str3, '\n')
  str
}
