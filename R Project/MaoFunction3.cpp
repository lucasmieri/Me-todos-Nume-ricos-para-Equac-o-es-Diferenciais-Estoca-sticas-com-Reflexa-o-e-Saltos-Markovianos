#include <Rcpp.h>
using namespace Rcpp;

double MaoSim(double Dt, unsigned int L, double X0,
              NumericVector mu, NumericVector sigma,
              IntegerVector A, NumericVector B)
{
    double X = X0;
    for(unsigned int j = 1; j < L+1; ++j)
    {
        X = X + Dt*mu[A[j]]*X + sigma[A[j]]*X*B[j-1];
        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
NumericVector MaoFuncRcpp(double Dt, unsigned int L, unsigned int M,
                          IntegerVector &A, NumericMatrix &B)
{
    //Global env data
    Environment env = Environment::global_env();
    double X0            = env["X0"];
    NumericVector mu     = env["mu"];
    NumericVector sigma  = env["sigma"];
   
    //Running Simulations 
    NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
        X[j] = MaoSim(Dt,L,X0,mu,sigma,A,B(j,_));
    
    return(X);
}

