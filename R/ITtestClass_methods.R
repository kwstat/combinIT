#' @export
print.ITtest <- function(x,...){
  msg1 <- paste(" Test:\t", x$test, '\n')
  msg2 <- paste("Data:\t", x$data.name, '\n')
  msg3 <- paste("Statistic = ", round(x$statistic, 3), '\n')
  msg4 <- paste("Exact Monte Carlo P-value = ", ifelse(is.numeric(x$pvalue.exact), round(x$pvalue.exact, 3), "NULL"), '\n')
  msg5 <- paste("Approximate P-value = ", ifelse(is.numeric(x$pvalue.appro), round(x$pvalue.appro, 3), "NULL"),  '\n')
  msg6 <- paste("Nsim = ", x$nsim, '\n')
  cat(msg1, msg2, msg3, msg4, msg5, msg6)
}