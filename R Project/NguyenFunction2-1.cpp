// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
#define R_INTERRUPT 1000
#define N_PRINT_STATS 10000

double NguyenSim(double Dt, unsigned int L, double X0, NumericVector &p0,
              NumericVector &mu, NumericVector &sigma, NumericVector &lambda,
              IntegerVector &E)
{
    double X = X0;
    double t = 0;
    unsigned int A = sample(E,1,true,p0)[0];
    
    NumericVector norms = rnorm(L);
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
        X = X + Dt*X*mu[A] + X*norms[j]*sqrt(Dt)*sigma[A]
              + 0.5*sigma[A]*sigma[A]*X*(Dt*norms[j]*norms[j] - Dt);

        if(j % R_INTERRUPT == 0) Rcpp::checkUserInterrupt();
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
    long double mean = 0, var = 0, X;
    double CI; //NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
    {
        //X[j] = NguyenSim(Dt,L,X0,p0,mu,sigma,lambda,E);
        X = NguyenSim(Dt,L,X0,p0,mu,sigma,lambda,E);
        
        mean = mean + (X-mean)/(j+1.0);
        if(j > 0)
            var = var + (mean-X)*(mean-X)/(j+1.0) - var/(j);
        
        if(j > 0 && j % N_PRINT_STATS == 0 )
        {
            CI = Rf_qnorm5(0.995,0,1,true,false)*var/sqrt(j+1);
            printf("Nguyen stats: %0.7Lf +- %0.5Lf (DONE %.0f%%)\n",mean,CI,((float)j+1)/M*100);
        }
    }
    return(mean);
}
