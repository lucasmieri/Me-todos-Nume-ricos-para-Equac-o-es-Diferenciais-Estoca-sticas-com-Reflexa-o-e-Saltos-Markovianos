#include <Rcpp.h>
using namespace Rcpp;

double NguyenSim(double Dt, unsigned int L, double X0, double A0,
              NumericVector mu, NumericVector sigma, NumericVector lambda)
{
    double X = X0;
    int    A = A0;
    double t = 0;
    double tau = rexp(1,lambda[A])[0];
    
    NumericVector norms = rnorm(L+1);
    for(unsigned int j = 0; j < L; ++j)
    {
        t += Dt;
        if(tau <= t)
        {
            if( A == 0 ) A = 1;
            else         A = 0;
            tau += rexp(1,lambda[A])[0];
        }
        X = X + Dt*X*mu[A] + X*norms[j]*sqrt(Dt)*sigma[A]
              + 0.5*sigma[A]*sigma[A]*X*(Dt*norms[j]*norms[j] - Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
NumericVector NguyenFuncRcpp(double Dt, unsigned int L, unsigned int M)
{
    //Global env data
    Environment env = Environment::global_env();
    double X0            = env["X0"];
    NumericVector p0     = env["p0"];
    NumericVector lambda = env["lambda"];
    NumericVector mu     = env["mu"];
    NumericVector sigma  = env["sigma"];
    IntegerVector E      = {0,1};
    int A0               = sample(E,1,false,p0)[0];
   
    //Running Simulations 
    NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
        X[j] = NguyenSim(Dt,L,X0,A0,mu,sigma,lambda);
    
    return(X);
}
