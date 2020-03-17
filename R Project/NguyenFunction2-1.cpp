// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

double NguyenSim(double Dt, unsigned int L, double X0, NumericVector &p0,
              NumericVector &mu, NumericVector &sigma, NumericVector &lambda,
              IntegerVector &E)
{
    double X = X0;
    double t = 0;
    unsigned int A = sample(E,1,true,p0)[0];
    
    //NumericVector norms = rnorm(L);
    double norm;
    double tau = rexp(1,lambda[A])[0];
    for(unsigned int j = 0; j < L; ++j)
    {
        t += Dt;
        if(tau <= t)
        {
            if( A == 0 ) A = 1;
            else         A = 0;
            tau += rexp(1,lambda[A])[0];
        }
        norm = rnorm(1)[0];
        X = X + Dt*X*mu[A] + X*norm*sqrt(Dt)*sigma[A]
              + 0.5*sigma[A]*sigma[A]*X*(Dt*norm*norm - Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
double NguyenFuncRcpp(double Dt, unsigned int L, unsigned int M)
{
    //Global env data
    Environment env = Environment::global_env();
    double X0            = env["X0"];
    NumericVector p0     = env["p0"];
    NumericVector lambda = env["lambda"];
    NumericVector mu     = env["mu"];
    NumericVector sigma  = env["sigma"];
    IntegerVector E      = {0,1};
   
    //Running Simulations 
    double mean = 0, alpha; //X(M);
    for(unsigned int j = 0; j < M; ++j)
    {
        alpha = ((double)j)/(j+1);
        mean = (alpha)*mean + (1-alpha)*NguyenSim(Dt,L,X0,p0,mu,sigma,lambda,E);
        //X[j] = NguyenSim(Dt,L,X0,p0,mu,sigma,lambda,E);
    }
    return(mean);
}
