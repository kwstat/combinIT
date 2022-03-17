#' @importFrom stats quantile
completed <- function(nsim) {
  out <- list(
    j = 1,
    pr = floor(quantile(1:nsim)),
    lab = c("0%", "25%", "50%", "75%", "100%"),
    nsim = nsim
  )
  structure(out, class = "completed")
}

nextc <- function(x, i) {
  UseMethod("nextc")
}

nextc.completed <- function(x, i) {
  message(x$lab[x$j], " completed")
  x$j <- x$j + 1
  return(x)
}
