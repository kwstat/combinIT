#' Boik's LBI test
#'
#' Computes LBI's test statistics \eqn{(tr(R'R))^2/(p tr((R'R)^2))}, where R is residual
#' matrix of input data matrix under additivity assumption, and returns corresponung p-value.
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects.
#' @param nsim Number of simulation for compueting exact p-value. The defaut value is 1000.
#' @return Exact and asymptotic p-value for input
#' @author Zahra. Shenavari, ....
#' @references Boik, R. J. (1993a). Testing additivity in two-way classifications
#'  with no replications: the locally best invariant test. Journal of Applied
#'  Statistics 20(1): 41-55.
#' @examples \dontrun{this is an example}
#' data(MVGH)
#' Boik.test(MVGH, nsim=10000)
#' @export
Boik.test <- function(x, nsim=1000, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)
    statistics <-Bfc(x,bl,tr,p)
    simu <- Bfsim(nsim,bl,tr,p)
    boik.p <- mean(statistics > simu)
    Tb<-(1/statistics-1)
    T<-p*q*Tb/2
    df<-(p+2)*(p-1)/2
    if(p==2) asyboik.p <-pbeta(Tb,1,(q-1)/2)
    else asyboik.p <- pchisq(T, df )
    out <- list(exact.pvalue = boik.p,asy.pvalue = asyboik.p,
                nsim = nsim,
                statistics = statistics)
  }
    return(out)
  }


#' @export
Boik.test.old <- function(x, nsim=1000, ...) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)

    statistics <-B.f(x,p)
    simu <-rep(0,0)
    for (i in 1:nsim){
      simu[i]<-B.f(matrix(rnorm(n),nrow = bl),p)
      #        cat(paste(round(i / nsim * 100), '% completed'))
      #        #Sys.sleep(.0001)
      #        if (i == nsim) cat(': Done')
      #        else cat('\014')
    }
    boik.p <- mean(statistics > simu)
    Tb<-(1/statistics-1)
    T<-p*q*Tb/2
    df<-(p+2)*(p-1)/2
    if(p==2) asyboik.p <-pbeta(Tb,1,(q-1)/2)
    else asyboik.p <- pchisq(T, df )
    out <- list(exact.pvalue = boik.p,asy.pvalue = asyboik.p,
                nsim = nsim,
                statistics = statistics)
    return(out)
  }
}

#' Malik's test for interaction
#'
#' Repors the exact p-value from \code{\link{Malik et al. (2016)}}.They proposed to partition
#' the residuals(components of residual matrix R of input data matrix under additivity
#' assumption) into three clusters using a suitable clustering method like “k-means clustering”.
#' The hypothesis of no interaction can be interpreted as the effects of the three
#' clusters are equal.Compute Malik's test statistics and corresponding p-value by Monte Carlo simulation.
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects
#' @param nsim Number of simulation for compueting exact p-value. The defaut value is 500
#' @return An exact p-value for input
#' @author Zahra. Shenavari, ..
#' @references Malik, W. A., Mo ̈hring, J., Piepho, H. P. (2016). A
#' clustering-based test for non-additivity in an unreplicated two-way layout.
#' Communications in Statistics-Simulation and Computation 45(2):660-670.
#' @examples \dontrun{this is an example}
#' data(impurity)
#' Malik.test(impurity,nsim=1000)
#' @export
Malik.test <- function(x, nsim=500) {
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    statistics <- M_f(x)
    simu <-rep(0,0)
    for (i in 1:nsim){
         y<-rnorm(n)
        simu[i] <- M_f(matrix(y,nrow = bl))
      }
      malik <- mean(statistics < simu)
    list(pvalue = malik,nsim=nsim,statistics=statistics)
  }
}

#' Kharrati-Kopaei and Miller's test for interaction
#'
#' Computes the test statistics and corresponding p-value based on inspecting all
#' Pairwise Interaction Contrasts (PICs) for testing interaction
#' proposed by Kharrati-Kopaei and Miller(2016).
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects
#' @param nsim Number of simulation for compueting exact p-value. The defaut value is 1000
#' @return An exact p-value for input
#' @author Zahra. Shenavari, ...
#' @references Kharrati-Kopaei, M., Miller, A. (2016). A method for testing interaction in
#'  unreplicated two-way tables: using all pairwise interaction contrasts. Statistical
#'   Computation and Simulation 86(6):1203-1215.
#' @examples \dontrun{this is an example}
#' data(ratio)
#' PIC.test(ratio,nsim=10000,nc0=10000)
#' @export
PIC.test <- function(x, nsim=1000,nc0=10000, ...) {

  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    kp<- kpr(bl, tr)
    c0<-C0(kp,n,nc0)
    statistics <-picf(y,kp,c0)
    simu <- PICfsim(nsim,kp,c0,n)
    PIC <- mean(statistics < simu)
    list(pvalue = PIC,nsim=nsim,statistics=statistics)
  }
}


#' @export
PIC.test.old <- function(x, nsim=1000,nc0=10000, ...) {

  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    n <- tr * bl
    kp<- kpr(bl, tr)
    c0<-mean(replicate(nc0,{median(abs(kp%*%rnorm(n)))}))
    statistics <-pic.f(y,kp,c0)
    simu <-rep(0,0)
    for (i in 1:nsim){
      simu[i]<-pic.f(rnorm(n),kp,c0)
#      cat(paste(round(i / nsim * 100), '% completed'))
      #Sys.sleep(.05)
#      if (i == nsim) cat(': Done')
#      else cat('\014')
    }
    PIC <- mean(statistics < simu)

    list(pvalue = PIC,nsim=nsim,statistics=statistics)
  }
}

#' Piepho's third version test for interaction
#'
#' Piepho (1994) proposed three test statistics.The third one is
#' based on Grubbs’ (1948) type estimator of variance for each level of block effect.
#' reports Piepho's test statistics and coresponding p-value by Monte Carlo datasets
#' or asympthotic distribution.
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects
#' @param nsim Number of simulation for compueting exact p-value. The defaut value is 1000
#' @return Exact and asymptotic p-value for input
#' @author Zahra. Shenavari, ...
#' @references Kharrati-Kopaei, M., Miller, A. (2016). A method for testing
#' interaction in unreplicated two-way tables: using all pairwise interaction contrasts.
#'  Statistical Computation and Simulation 86(6):1203-1215.
#' @examples \dontrun{this is an example}
#' data(MVGH)
#' piepho.test(MVGH,sim=1000)
#' @export
piepho.test<-function(x,nsim=1000,...){
  if(!is.matrix(x)){
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n<-tr * bl
    statistics <-piephoC(x,bl,tr)
    simu <- Piephosim(nsim,bl,tr)
    pieph <- mean(statistics < simu)
    df = bl-1
    asypieph <- 1-pchisq(statistics, df = df)
    out <- list(exact.pvalue = pieph,asypvalue = asypieph,
                  nsim = nsim
                  ,statistics = statistics)
    return(out)
  }
}


#' @export
piepho.test.old<-function(x,nsim=1000,...){
  if(!is.matrix(x)){
    stop("The input should be a matrix")
  } else {
    tr <- ncol(x)
    bl <- nrow(x)
    n<-tr * bl
    statistics <-piepho(x,bl,tr)
    simu <-rep(0,0)
    for (i in 1:nsim){
      simu[i]<-piepho(matrix(rnorm(n),nrow = bl),bl,tr)
      #cat(paste(round(i / nsim * 100), '% completed'))
      #Sys.sleep(.05)
      #if (i == nsim) cat(': Done')
      #else cat('\014')
    }
    pieph <- mean(statistics < simu)
    df = bl-1
    asypieph <- 1-pchisq(statistics, df = df)
    
    out <- list(exact.pvalue = pieph,asypvalue = asypieph,
                nsim = nsim
                ,statistics = statistics)
    return(out)
  }
}

#' Kharrati-Kopaei and Sadooghi-Alvandi's test for interaction
#'
#' computes Kharrati-Kopaei and Sadooghi-Alvandi's test statistics and corresponding p-value.
#' They proposed a test procedure for testing non-additivity based on splitting the b-by-t data
#' matrix into two non empty sub-matrixes by dividing the rows, each having at least two rows.
#' Suppose that b≥t and b≥4.Divided mean square errors for any two sub-matrixes such that
#' maximize its ratio and computed the corresponding p-value. Their test statistics are the
#' minimum value of p-values over all possible configurations: 2^(b-1)-b-1.
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects
#' @details If rows numer(b) of data matrix is less than it's columns number(t) at fiest
#'  transposing data matrix. In addition requires that data matrix has more than three
#'   rows or columns
#' @param nsim Number of simulation for compueting exact p-value
#' @param dist If dist="sim", we used Monte Carlo simulation for estimating exact p-value,
#'  and if dist="asy", Bonferroni-adjusted p-value computed. The defaut value is "sim"
#' @return A p-value for input
#' @author Zahra. Shenavari, ...
#' @references Kharrati-Kopaei, M., Sadooghi-Alvandi, S. M. (2007). A New Method for
#'  Testing Interaction in Unreplicated Two-Way Analysis of Variance. Communications
#'   in Statistics-Theory and Methods 36:2787–2803.
#' @examples \dontrun{this is an example}
#' data(impurity)
#' KKSA.test(impurity,nsim=1000,dist = "sim")
#' @export
KKSA.test<-function(x,nsim=1000,distr = "sim",...){

  if(!is.matrix(x)){
    stop("The input should be a matrix")
  }else{
    bl <- nrow(x)
    tr <- ncol(x)
    n<-tr * bl
    if (bl < tr)
    {warning("transpose the input matrix")
      x<-t(x);te<-bl;bl<-tr;tr<-te}
    if (bl < 4) {
      stop("KKSA needs at least 4 levels of blocking factor")

    } else{
      cck <- 2^(bl - 1) - 1 - bl
      statistics <-kk_f(x)
      if(distr != "sim" && distr != "asy") distr="sim"

      if (distr == "sim")
       {
        simu <-rep(0,0)
        for (i in 1:nsim){
          simu[i]<-kk_f(matrix(rnorm(n),nrow=bl))
        }
        KKSA.p <- mean(statistics > simu)
      } else if (distr == "asy") {
        KKSA.p<-statistics*cck
        KKSA.p<- min(1,KKSA.p)
       }
      out <- list(pvalue = KKSA.p,
                nsim = nsim,distr=distr,
                statistics = statistics)
      return(out)
    }
  }
}

#' Franck et al.'s test for interaction
#'
#' computes Franck et al's. (2013) test statistic,ACMIF, and corresponding p-value.
#' @param x A b-by-t data matrix, which rows corresponding
#'  to b-block effects and columns are t-treatment effects
#' @details If rows numer(b) of data matrix is less than it's columns number(t) at fiest
#'  transposing data matrix. In addition requires that data matrix has more than two
#'   rows or columns
#' @param nsim Number of simulation for compueting exact p-value
#' @param dist If dist="sim", we used Monte Carlo simulation for estimating exact p-value,
#'  and if dist="asy", Bonferroni-adjusted p-value computed. The defaut value is "sim"
#' @return A p-value for input
#' @author Zahra. Shenavari, ...
#' @references Franck, C., Nielsen, D., Osborne, J. A. (2013). A method for detecting
#' hidden additivity in two-factor unreplicated experiments. Computational Statistics
#'  and Data Analysis 67:95-104.
#' @examples \dontrun{this is an example}
#' data(cnv6)
#' hiddenf.test(cnv6,nsim=1000,dist = "sim")
#' @export
hiddenf.test<-function(x,nsim=1000,dist = "sim",...){

  if(!is.matrix(x)){
    stop("The input should be a matrix")
  } else {
    bl <- nrow(x)
    tr <- ncol(x)
    n<-tr * bl
    if (bl < tr) {warning("transpose the input matrix")
      x<-t(x);te<-bl;bl<-tr;tr<-te}
    if (bl < 3) {
      stop("hiddenf needs at least 3 levels of blocking factor")

    }else{
     cch <- 2^(bl - 1) - 1
     statistics <-hh.f(x,bl)
     if(dist != "sim" || dist != "asy") dist="sim"

     if (dist == "sim")
     {
       simu <-rep(0,0)
       for (i in 1:nsim){
         simu[i]<-hh.f(matrix(rnorm(n),nrow=bl),bl)
         cat(paste(round(i / nsim * 100), '% completed'))
         #Sys.sleep(.1)
         if (i == nsim) cat(': Done')
         else cat('\014')
       }
       hidden <- mean(statistics < simu)
     } else if (dist == "asy") {
       adjpvalue<-(1-pf(statistics,(tr-1),(tr-1)*(bl-2)))*cch
       hidden<- min(1,adjpvalue)
     }
     out <- list(pvalue = hidden,
                nsim = nsim,dist=dist
                ,statistics = statistics)
     return(out)
     }
   }
}

#' Combined several interaction tests
#'
#' Reports p-values tests for non-additivity developed by Boik (1993a), Piepho (1994),
#' Kharrati-Kopaei and Sadooghi-Alvandi (2007), Franck et al. (2013), Malik et al. (2016)
#' and Kharrati-Kopaei and Miller (2016). In addition by use of four combination methods:
#' Bonferroni, Sidak, Jacobi expantion, and Gaussian Copula combined reported p-values.
#' @param x A data matrix in two-factor analysis
#' @details If rows numer(b) of data matrix is less than it's columns number(t) we
#'  transpose data matrix. In addition requires that data matrix has more than two
#'   rows or columns.
#' @param nsim Number of simulation for compueting exact p.value
#' @param dist If dist="sim", we used Monte Carlo simulation for estimating exact p-value,
#'  and if dist="asy", the p.values is
#' estimated from an asymptotic distribuion. The defaut value is "sim".
#' @return  A p.value for input
#' @details Needs "mvtnorm" packages
#' @author Zahra. Shenavari, ...
#' @references Shenavari, Z.,Kharrati-Kopaei, M. (2018). A Method for Testing Additivity
#' in Unreplicated Two-Way Layouts Based on Combining Multiple Interaction Tests.International
#' Statistical Review
#' @examples \dontrun{this is an example}
#' data(cnv6)
#' combinep(cnv6,nsim=500,nc0=10000)
#' @export
combinep<-function(x,nsim=500,nc0=10000,...){
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    y <- c(t(x))
    tr <- ncol(x)
    bl <- nrow(x)
    if (bl < tr) warning("transpose the input matrix")
    x<-t(x);te<-bl;bl<-tr;tr<-te}
    if (bl < 3) {
      stop("hiddenf needs at least 3 levels of blocking factor")

    }else{
    n <- tr * bl
    block <- gl(bl, tr)
    treatment <- gl(tr, 1, bl * tr)
    p <- min(tr - 1, bl - 1)
    q <- max(tr - 1, bl - 1)
    cck <- 2^(bl - 1) - 1 - bl
    cch <- 2^(bl - 1) - 1
    kp <- kpr(bl, tr)
    c0<-mean(replicate(nc0,{median(abs(kp%*%rnorm(n)))}))

    sta<-bmp.f(x,y, block , treatment,bl,tr,p)

    Bstat <-sta$Boik
    Mstat <-sta$Tc
    pistat<-sta$piepho
    pstat <-pic.f(y,kp,c0)
    if(bl==3)Hstat<-hh.f(x,bl)
    else{
      Ksimu <-rep(0,0)
      kh<-kh.f(x,bl,tr)
      Kstat<-kh$fmin
      Hstat<-kh$fmax
    }

      Bsimu <-Msimu <-psimu <-pisimu <-Hsimu <-rep(0,0)
      for (i in 1:nsim){
        y<-rnorm(n)
        x<-matrix(y,nrow = bl,byrow=TRUE)
        sta<-bmp.f(x,y, block , treatment,bl,tr,p)
        Bsimu[i]<-sta$Boik
        Msimu[i]<-sta$Tc
        pisimu[i]<-sta$piepho
        psimu[i]<-pic.f(y,kp,c0)
        if(bl==3)
          Hsimu[i]<-hh.f(x,bl)else{
        kh<-kh.f(x,bl,tr)
        Ksimu[i]<-kh$fmin
        Hsimu[i]<-kh$fmax
          }
        cat(paste(round(i / nsim * 100), '% completed'))
        Sys.sleep(.1)
        if (i == nsim) cat(': Done')
        else cat('\014')
      }
      Boik.pvalue <- mean(Bstat > Bsimu)
      piepho.pvalue <- mean(pistat < pisimu)
      PIC.pvalue <- mean(pstat < psimu)
      Malik.pvalue <- mean(Mstat < Msimu)
      hiddenf.pvalue <- mean(Hstat < Hsimu)
      if(bl==3)
        KKSA.pvalue<- NULL else{

        KKSA.pvalue<- mean(Kstat > Ksimu)
      }
    pvalues <- c(Boik.pvalue,piepho.pvalue,hiddenf.pvalue,Malik.pvalue,PIC.pvalue,KKSA.pvalue)
    cp<-comb(pvalues)
    Bonferroni<-cp$Bon
    GC<-cp$GC
    Sidak<-cp$Sidak
    jacobi<-cp$jacobi
    list(nsim=nsim,piepho.pvalue=piepho.pvalue,Boik.pvalue=Boik.pvalue
         ,Malik.pvalue=Malik.pvalue,PIC.pvalue=PIC.pvalue
         ,KKSA.pvalue=KKSA.pvalue,hiddenf.pvalue=hiddenf.pvalue,
         Bonferroni=Bonferroni,Sidak=Sidak,jacobi=jacobi,GC=GC)
    }
  }


#' Interaction plot
#'
#' @param x A data matrix in two-factor analysis
#' @return  An interaction plot for input
#' @author Zahra. Shenavari, ...
#' @examples \dontrun{this is an example}
#' data(cnv6)
#' interactionplot(cnv6)
#' @export
interactionplot<-function(x,...){
  if (!is.matrix(x)) {
    stop("The input should be a matrix")
  } else {
    par(mfcol=c(1,2))
    t <- ncol(x)
    b <- nrow(x)
  matplot(t(x), type = "l",xaxt = "n", ylab = "y", xlab = "Factor1(column)",lty = 1:b,...)
  axis(1, at = 1:t, labels = 1:t, cex.axis = 1)
  matplot(x, type = "l",xaxt = "n", ylab = "", xlab = "Factor2(row)",lty = 1:t,...)
  axis(1, at = 1:b, labels = 1:b, cex.axis = 1)
  }
}














