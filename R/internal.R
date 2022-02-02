

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

#' internal function for Malik method
#'
#' @keywords internal
#'
M.f <- function(x, y, block, treatment) {
  RES <- t(t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))
  r <- c(t(RES))
  kmean <- kmeans(x = r, centers = 3, nstart = 100)
  af <- kmean$cluster
  modclus <- lm(y ~ block + treatment + as.factor(af))
  amodclus <- anova(modclus)
  Tc <- (amodclus[3, 2] / amodclus[3, 1]) / (amodclus[4, 2] / amodclus[4, 1])
  return(Tc)
}

#' This is an internal function which compute Kronecker product for PIC method
#'
#' @keywords internal
#'
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

#' This is an internal function for PIC method
#'
#' @keywords internal
#'
pic.f <- function(y, kp, c0) {
  z <- abs(kp %*% y)
  s0 <- median(z) / c0
  PSE <- median(z[z <= 5 * s0])
  PIC <- max(z) / PSE
  return(PIC)
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

#' internal function for KKSA method
#'
#' @param x A data matrix
#' @keywords internal
#'
kk.f <- function(x, bl, tr) {
  Nrow <- 2:(as.integer(bl / 2))
  fvalues <- rep(0, 0)
  pvalues <- rep(0, 0)
  count <- 0
  for (i in Nrow) {
    ind <- combn(bl, i)
    Nsplit <- ncol(ind)
    if (bl / 2 == i) Nsplit <- Nsplit / 2
    for (j in 1:Nsplit) {
      count <- count + 1
      yb1 <- x[ind[, j], ]
      yb2 <- x[-c(ind[, j]), ]
      rss1 <- sum((t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
      rss2 <- sum((t(yb2 - apply(yb2, 1, mean) + mean(yb2)) - apply(yb2, 2, mean))^2)
      dfn <- (tr - 1) * (i - 1)
      dfd <- (bl - i - 1) * (tr - 1)
      fvalues[count] <- (rss1 * (bl - i - 1)) / (rss2 * (i - 1))
      if (fvalues[count] < 1) fvalues[count] <- 1 / fvalues[count]
      pvalues[count] <- 1 - pf(fvalues[count], dfn, dfd) + pf(1 / fvalues[count], dfn, dfd)
    }
  }
  KKSA <- min(pvalues)
  return(KKSA)
}

#' internal function for hidden method
#'
#' @param x A data matrix
#' @keywords internal
#'
hh.f <- function(x, bl) {
  Nrow <- 2:(as.integer(bl / 2))
  sse <- sum((t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))^2)
  hvalues <- rep(0, 0)
  count <- 0
  for (i in Nrow) {
    ind <- combn(bl, i)
    Nsplit <- ncol(ind)
    if (bl / 2 == i) Nsplit <- Nsplit / 2
    for (j in 1:Nsplit) {
      count <- count + 1
      yb1 <- x[ind[, j], ]
      yb2 <- x[-c(ind[, j]), ]
      rss1 <- sum((t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
      rss2 <- sum((t(yb2 - apply(yb2, 1, mean) + mean(yb2)) - apply(yb2, 2, mean))^2)
      sse7 <- rss1 + rss2
      hvalues[count] <- (sse - sse7) * (bl - 2) / sse7
    }
  }
  for (d in 1:bl) {
    count <- count + 1
    yb1 <- x[-d, ]
    sse7 <- sum((t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
    hvalues[count] <- (sse - sse7) * (bl - 2) / sse7
  }

  fmax <- max(hvalues)
  return(fmax)
}

#' internal function for hidden and KKSA methods
#'
#' @param x A data matrix
#' @keywords internal
#'
kh.f <- function(x, bl, tr) {
  Nrow <- 2:(as.integer(bl / 2))
  sse <- sum((t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))^2)
  hvalues <- rep(0, 0)
  fvalues <- rep(0, 0)
  pvalues <- rep(0, 0)
  count <- 0
  for (i in Nrow) {
    ind <- combn(bl, i)
    Nsplit <- ncol(ind)
    if (bl / 2 == i) Nsplit <- Nsplit / 2
    for (j in 1:Nsplit) {
      count <- count + 1
      yb1 <- x[ind[, j], ]
      yb2 <- x[-c(ind[, j]), ]
      rss1 <- sum((t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
      rss2 <- sum((t(yb2 - apply(yb2, 1, mean) + mean(yb2)) - apply(yb2, 2, mean))^2)
      sse7 <- rss1 + rss2
      hvalues[count] <- (sse - sse7) * (bl - 2) / sse7
      #--------------------
      dfn <- (tr - 1) * (i - 1)
      dfd <- (bl - i - 1) * (tr - 1)
      fvalues[count] <- (rss1 * (bl - i - 1)) / (rss2 * (i - 1))
      if (fvalues[count] < 1) fvalues[count] <- 1 / fvalues[count]
      pvalues[count] <- 1 - pf(fvalues[count], dfn, dfd) + pf(1 / fvalues[count], dfn, dfd)
    }
  }
  fmin <- min(pvalues)
  for (d in 1:bl) {
    count <- count + 1
    yb1 <- x[-d, ]
    sse7 <- sum((t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
    hvalues[count] <- (sse - sse7) * (bl - 2) / sse7
  }

  fmax <- max(hvalues)
  list(fmin = fmin, fmax = fmax)
}

#' internal function fo rBoik.Malik.Piepho
#'
#' @keywords internal
#'
bmp.f <- function(x, y, block, treatment, bl, tr, p) {
  RES <- t(t(x - apply(x, 1, mean) + mean(x)) - apply(x, 2, mean))
  r <- c(t(RES))
  #---------------
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
  #---------------
  kmean <- kmeans(x = r, centers = 3, nstart = 100)
  assn <- kmean$cluster
  modclus <- lm(y ~ block + treatment + as.factor(assn))
  amodclus <- anova(modclus)
  Tc <- (amodclus[3, 2] / amodclus[3, 1]) / (amodclus[4, 2] / amodclus[4, 1])
  #-----------------
  EE1 <- crossprod(RES)
  EE2 <- EE1 %*% EE1
  trace1 <- sum(diag(EE1))
  trace2 <- sum(diag(EE2))
  Boik <- trace1^2 / (p * trace2)
  #-----------
  list(Boik = Boik, Tc = Tc, piepho = piepho)
}

#' internal function forcombining pvalues
#' #library(mvtnorm)
#'
#' @keywords internal
#'
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
