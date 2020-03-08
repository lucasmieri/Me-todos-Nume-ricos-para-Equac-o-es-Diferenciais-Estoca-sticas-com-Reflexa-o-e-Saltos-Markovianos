#include <Rcpp.h>
using namespace Rcpp;

double MaoSim(double Dt, unsigned int L, double X0, double A0,
              NumericVector mu, NumericVector sigma,
              NumericMatrix P, IntegerVector E)
{
    double X = X0;
    int    A = A0;
   
    NumericVector norms = rnorm(L);
    for(unsigned int j = 0; j < L; ++j)
    {
        A = sample(E,1,false,(NumericVector)P(A,_))[0];
        X = X + Dt*mu[A]*X + sigma[A]*X*norms[j]*sqrt(Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
NumericVector MaoFuncRcpp(double Dt, unsigned int L, unsigned int M)
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
   
    //Gerar P (usar expm)
    NumericMatrix Rho(2,2);
    Rho(0,0) = -lambda[0]; Rho(0,1) = lambda[0];
    Rho(1,0) = lambda[1];  Rho(1,1) = -lambda[1];
    Function expm("expm"); 
    NumericMatrix P = expm(Rho*Dt);
   
    Function set_seed("set.seed");
    
    //Running Simulations 
    NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
    {
        set_seed(j); 
        X[j] = MaoSim(Dt,L,X0,A0,mu,sigma,P,E);
    }
    
    return(X);
}

