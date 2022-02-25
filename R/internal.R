

#' This is an internal function for Boik method
#'
#' @keywords internal
#'
B.f <- function(x, p) {
  RES <- t(t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))
  EE1 <- crossprod(RES)
  EE2 <- EE1 %*% EE1
  trace1 <- sum(diag(EE1))
  trace2 <- sum(diag(EE2))
  Boik <- trace1^2 / (p * trace2)
  return(Boik)
}


#' This is an internal function which compute Kronecker product for PIC method
#'
#' @keywords internal
#'
#'@importFrom utils combn
kpr <- function(bl, tr) {
  wa <- combn(bl, 2)
  wb <- combn(tr, 2)
  cb <- matrix(0, nrow = choose(bl, 2), ncol = bl)
  ct <- matrix(0, nrow = choose(tr, 2), ncol = tr)
  for (i in 1:choose(bl, 2)) {
    cb[i, wa[1, i]] <- 1
    cb[i, wa[2, i]] <- -1
  }
  for (i in 1:choose(tr, 2)) {
    ct[i, wb[1, i]] <- 1
    ct[i, wb[2, i]] <- -1
  }
  c <- kronecker(cb, ct)
  return(c)
}


#' internal function for Piepho method
#'
#' @keywords internal
#'
piepho <- function(x, bl, tr) {
  RES <- t(t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))
  W <- apply(RES^2, 1, sum)
  delta <- (bl * (bl - 1) * W - sum(W))
  h1 <- 0
  for (i in 1:(bl - 1)) {
    for (j in (i + 1):bl) {
      h1 <- (delta[i] * delta[j]) + h1
    }
  }
  U <- 2 * bl * h1 / ((bl - 1) * (sum(delta)^2))
  piepho <- -(tr - 1) * (bl - 1) * (bl - 2) * log(U) / 2
  return(piepho)
}



#' internal function forcombining pvalues

#'
#' @keywords internal
#'
#' @importFrom mvtnorm pmvnorm

comb <- function(pvalues) {
  P <- pvalues
  P[P == 0] <- 10^(-6)
  P[P == 1] <- 1 - 10^(-6)
  k <- length(P)
  T <- qnorm(P)
  r <- 1 - var(T)
  rohat <- max(-1 / (k - 1), r)
  j <- matrix(1, nrow = k, ncol = k)
  i <- diag(k)
  S <- (1 - rohat) * i + rohat * j
  minp <- min(P)
  m <- qnorm(minp)
  GC <- 1 - pmvnorm(lower = rep(m, k), upper = Inf, sigma = S)
  q0 <- max(1 / minp - 1, 1)
  q <- min(q0, (k - 1))
  Bon <- min(minp * k, 1)
  jacobi <- 1 - (1 - minp)^q # minpv~Betha(1,q)
  Sidak <- 1 - (1 - minp)^k # minpv~Betha(1,k)
  list(Bon = Bon, Sidak = Sidak, jacobi = jacobi, GC = GC) # GC=Goussian Copula=MVNormal
}
