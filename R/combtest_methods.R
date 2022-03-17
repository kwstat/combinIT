#' @export
print.combtest <- function(x, ...) {
  msg1 <- paste(" Test:\t", x$test, "\n")
  msg2 <- paste("Data:\t", x$data.name, "\n")
  msg3 <- paste("Piepho Test: Statistics = ", round(x$Piepho.Stat, 5), ", Pvalue = ", round(x$Boik.pvalue, 5),"\n")
  msg4 <- paste("Boik Test: Statistics = ", round(x$Boik.Stat, 5),", Pvalue = ", round(x$Boik.pvalue, 5),"\n")
  msg5 <- paste("Malik Test: Statistics = ", round(x$Malik.Stat, 5),", Pvalue = ", round(x$Malik.pvalue, 5),"\n")
  msg6 <- paste("KKM Test: Statistics = ", round(x$KKM.Stat, 5),", Pvalue = ", round(x$KKM.pvalue, 5),"\n")
  msg7 <- paste("KKSA Test: Statistics = ", round(x$KKSA.Stat, 5),", Pvalue = ", round(x$KKSA.pvalue, 5),"\n")
  msg8 <- paste("Malik Test: Statistics = ", round(x$Franck.Stat, 5),", Pvalue = ", round(x$Franck.pvalue, 5),"\n")
  msg9 <- paste("Bonferroni method: Pvalue =", round(x$Bonferroni, 5),"\n")
  msg10 <- paste("Sidak method: Pvalue =", round(x$Sidak, 5),"\n")
  msg11 <- paste("Jacobi method: Pvalue =", round(x$Jacobi, 5),"\n")
  msg12 <- paste("Gaussian copula: Pvalue =", round(x$GC, 5),"\n")
  cat(msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, msg10, msg11, msg12)
}
