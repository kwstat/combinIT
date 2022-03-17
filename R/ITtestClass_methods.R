#' @export
print.ITtest <- function(x, ...) {
  msg1 <- paste(" Test:\t", x$test, "\n")
  msg2 <- paste("Data:\t", x$data.name, "\n")
  msg3 <- paste("Statistics = ", round(x$statistic, 5), "\n")
  msg4 <- paste("P-value = ", x$pvalue, "\n")
  msg5 <- paste("Nsim = ", x$nsim, "\n")
  msg6 <- paste("dist = ", x$dist, "\n")
  cat(msg1, msg2, msg3, msg4, msg5, msg6)
}
