// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <random>
#include <cmath>
using namespace Rcpp;

double MaoSim(double Dt, unsigned int L, double X0, NumericVector &p0,
              NumericVector &mu, NumericVector &sigma,
              NumericMatrix &P, IntegerVector &E)
{
    double X = X0;
    unsigned int A = sample(E,1,true,p0)[0];
   
    //NumericVector norms = rnorm(L);
    for(unsigned int j = 0; j < L; ++j)
    {
        A = sample(E,1,true,(NumericVector)P(A,_))[0];
        X = X + Dt*mu[A]*X + sigma[A]*X*rnorm(1)[0]*sqrt(Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
double MaoFuncRcpp(double Dt, unsigned int L, unsigned int M)
{
    //Global env data
    Environment env = Environment::global_env();
    double X0            = env["X0"];
    NumericVector p0     = env["p0"];
    NumericVector lambda = env["lambda"];
    NumericVector mu     = env["mu"];
    NumericVector sigma  = env["sigma"];
    IntegerVector E      = {0,1};
   
    //Gerar P (usar expm)
    NumericMatrix Rho(2,2);
    Rho(0,0) = -lambda[0]; Rho(0,1) = lambda[0];
    Rho(1,0) = lambda[1];  Rho(1,1) = -lambda[1];
    Function expm("expm"); 
    NumericMatrix P = expm(Rho*Dt);
   
    //Running Simulations 
    double mean = 0, alpha; //X(M);
    for(unsigned int j = 0; j < M; ++j)
    {
        alpha = ((double)j)/(j+1);
        mean = (alpha)*mean + (1-alpha)*MaoSim(Dt,L,X0,p0,mu,sigma,P,E);
        //X[j] = MaoSim(Dt,L,X0,p0,mu,sigma,P,E);
    }
    
    return(mean);
}


// [[Rcpp::export]]
IntegerVector sample_cpp(int n, NumericVector p){
    IntegerVector E = {0,1};
    return sample(E,n,true,p);
}