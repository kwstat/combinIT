#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
//' @importFrom Rcpp sourceCpp
//' @useDynLib combinIT
//' 
// [[Rcpp::export]]
float Bfc(arma::mat x,int bl, int tr,int p) {
  arma::vec RowMean = arma::mean(x,1);
  arma::vec ColMean= trans(mean(x,0));
  double Mean = accu(x)/(tr*bl);
  arma::mat y(bl,tr),yt1(bl,tr),yt2(bl,tr);
  for(int i=0; i<bl;i++)
    {
      for(int j=0;j< tr;j++)
      {
       y(i,j) = x(i,j)-RowMean(i)-ColMean(j)+Mean; // Zahra: please correct it by yourself
      }
   }
 yt1 = y.t() * y;
 yt2 = yt1*yt1;
  float Boik = trace(yt1)*trace(yt1) / (p*trace(yt2));
  return Boik;
}
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::vec Bfsim(int nsim,int bl, int tr,int p){
  arma::mat sam(bl,tr);
  arma::vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(bl,tr);
    out(i)=Bfc(sam,bl,tr,p);
  }
  return out;
}
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
double picf(arma::vec y,arma::mat kp,float c0){
  arma::vec z= kp * y;
  for(unsigned int i=0;i<kp.n_rows;i++)
    z(i)=fabs(z(i));
    arma::vec s0=median(z,0);
    arma::uvec ids = find(z <= (5*s0(0)) );
    arma::vec PSE=median(z.elem(ids),0);
    arma::vec PIC=max(z,0)/PSE(0);
    return PIC(0);
}
//' Module  function
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::vec PICfsim(int nsim,arma::mat kp, float c0, int n){
  arma::vec sam(n);
  arma::vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(n);
    out(i)=picf(sam,kp,c0);
  }
  return out;
}
//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
double C0(arma::mat kp, int n,int nc0){
  arma::vec sim(nc0);
  arma::vec norsam(n);
  for(int i=0;i<nc0;i++)
  {
    arma::vec temp=kp*norsam.randn(n);
    for(unsigned int j=0;j<temp.n_rows;j++)
      temp(j) = fabs(temp(j));
    arma::vec me=median(temp,0);
    sim(i)=me(0);
  }
  arma::vec out = mean(sim,0);
  return out(0);
}
//' @importFrom Rcpp sourceCpp
//' @useDynLib combinIT
//' 
// [[Rcpp::export]]
double piephoC(arma::mat x,int bl, int tr) {
  arma::vec RowMean = arma::mean(x,1);
  arma::vec ColMean= arma::trans(arma::mean(x,0));
  double Mean = arma::as_scalar(arma::accu(x)/(tr*bl));
  arma::mat Res2(bl,tr);
  for(int i=0; i<bl;i++)
  {
    for(int j=0;j< tr;j++)
    {
      double Res = arma::as_scalar(x(i,j)-RowMean(i)-ColMean(j)+Mean);
      Res2(i,j) = Res*Res; 
    }
  }
  arma::vec RowSum = arma::sum(Res2,1);
  arma::vec delta = (bl*(bl-1)*RowSum-sum(RowSum));
  double h1 = 0;
  for(int i=0;i< (bl-1);i++)
    for(int j=i+1; j < bl ; j++)
      h1 += arma::as_scalar(delta(i)*delta(j));
    double U = arma::as_scalar(2*bl*h1/((bl-1)*pow(sum(delta),2)));
    double piepho = -(tr-1)*(bl-1)*(bl-2)*log(U)/2;
    return piepho;
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::vec Piephosim(int nsim,int bl, int tr){
  arma::mat sam(bl,tr);
  arma::vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(bl,tr);
    out(i)=piephoC(sam,bl,tr);
  }
  return out;
}

using namespace arma;
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
double  M_f(arma::mat x) {            // Hossein's Codes....
  int bl = x.n_rows;
  int tr = x.n_cols;
  int n = bl*tr;
  arma::vec treatment = arma::repelem(arma::regspace(1,  bl), tr, 1);
  arma::vec y = arma::trans(x.as_row());
  arma::mat RES(bl,tr);
  arma::vec RowMean = arma::mean(x,1);
  arma::vec ColMean= trans(mean(x,0));
  double Mean = accu(x)/(tr*bl);
  for(int i=0; i<bl;i++)
  {
    for(int j=0;j< tr;j++)
    {
      RES(i,j) = x(i,j)-RowMean(i)-ColMean(j)+Mean; 
    }
  }
  arma::mat r = RES.as_row();  
  arma::mat centers;
  arma::kmeans(centers, r, 3, static_spread, 30, false);
  arma::vec af(n,fill::ones);
  for(int i=0; i<n;i++)
  {
    if(abs(r(i)-centers(0,1)) < abs(r(i)-centers(0,0)) && abs(r(i)-centers(0,1)) < abs(r(i)-centers(0,2)))
    {
      af(i)=2;
    } else if(abs(r(i)-centers(0,2)) < abs(r(i)-centers(0,0)) && abs(r(i)-centers(0,2)) < abs(r(i)-centers(0,1)))
    {
      af(i)=3;
    }
  }
  
  arma::vec  xi(n*3, fill::zeros), a1(n, fill::ones),a2(bl, fill::ones),a3(tr, fill::ones);
  arma::mat B3 = arma::diagmat(a3), B2 = arma::diagmat(a2);
  int j;
  for(int i=0; i<n;i++)
  {
    j = af(i);
    xi(3*i+j-1)= 1;
  }
  arma::mat  Xi = arma::trans( arma::reshape(xi, 3, n));
  arma::mat  K1 = arma::kron(a2,B3);
  arma::mat  K2 = arma::kron(B2,a3);
  arma::mat X = arma::join_horiz(a1,K1,K2,Xi);
  arma::vec yhat = X*arma::pinv(arma::trans(X)*X)*arma::trans(X)*y;
  double  SSE = arma::sum(arma::square(y-yhat));
  double Tc = ((arma::sum(arma::square(arma::vectorise(r)))-SSE)/2)/(SSE/((tr-1)*(bl-1)-2));
  return Tc;
}