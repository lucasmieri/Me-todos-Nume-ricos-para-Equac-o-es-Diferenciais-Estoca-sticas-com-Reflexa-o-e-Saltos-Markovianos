#include <Rcpp.h>
using namespace Rcpp;

double NguyenSim(double Dt, unsigned int L, double X0,
              NumericVector mu, NumericVector sigma,
              IntegerVector A, NumericVector B)
{
    double X = X0;
    for(unsigned int j = 1; j < L+1; ++j)
    {
        X = X + Dt*X*mu[A[j]] + X*B[j-1]*sigma[A[j]]
              + 0.5*sigma[A[j]]*sigma[A[j]]*X*(B[j-1]*B[j-1] - Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
NumericVector NguyenFuncRcpp(double Dt, unsigned int L, unsigned int M,
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
        X[j] = NguyenSim(Dt,L,X0,mu,sigma,A,B(j,_));
    
    return(X);
}
