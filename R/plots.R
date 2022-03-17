#' Interaction plot
#'
#' @param x numeric matrix, \eqn{b \times a} data matrix where the number of rows and columns are corresponding to the block and treatment levels, respectively.
#' @param ... plot parameters
#' @return  An interaction plot for input
#' @author Shenavari, Z.; Haghbin, H.; Kharrati-Kopaei, M.; Najibi, S.M.
#' @examples \dontrun{this is an example}
#' data(CNV)
#' interactionplot(CNV)
#' @importFrom graphics axis legend matplot par
#' @export
interactionplot <- function(x, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    t <- ncol(x)
    b <- nrow(x)
    oldpar <- par(mfcol=c(1,1))
    on.exit(par(oldpar))
    par(mfcol = c(1, 2))
    matplot(t(x), type = "b", xaxt = "n", ylab = "Observed values", xlab = "Factor1(column)", lty = 1:b, ...)
    axis(1, at = 1:t, labels = 1:t, cex.axis = 1)
    legend("topright", rep(paste0("row", 1:b)), lty = 1:b, bty = "n", cex = 0.7)
    matplot(x, type = "b", xaxt = "n", ylab = "", xlab = "Factor2(row)", lty = 1:t, ...)
    legend("topright", rep(paste0("col", 1:t)), lty = 1:t, bty = "n", cex = 0.7)
    axis(1, at = 1:b, labels = 1:b, cex.axis = 1)
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("combinIT", libpath)
}
