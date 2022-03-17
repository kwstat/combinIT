#' @export
print.comtest <- function(x,...){
  msg1 <- paste(" test:\t",x$test,'\n')
  msg2 <- paste("data:\t",x$data.name,'\n')
  msg3 <- paste("Statistics = ",round(x$statistic,2),'\n')
  msg4 <- paste("p-value = ",x$pvalue,'\n')
  msg5 <- paste("nsim = ",x$nsim,'\n')
  msg6 <- paste("dist = ",x$dist,'\n')
  cat(msg1,msg2,msg3,msg4,msg5,msg6)
}