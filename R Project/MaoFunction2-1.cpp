// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
#define R_INTERRUPT 1000
#define N_PRINT_STATS 10000

double MaoSim(double Dt, unsigned int L, double X0, NumericVector &p0,
              NumericVector &mu, NumericVector &sigma,
              NumericMatrix &P, IntegerVector &E)
{
    double X = X0;
    unsigned int A = sample(E,1,true,p0)[0];
   
    NumericVector norms = rnorm(L);
    NumericVector unifs = runif(L);
    for(unsigned int j = 0; j < L; ++j)
    {
        if(unifs[j] <= P(A,0)) A = 0;
        else                   A = 1;
        //A = sample(E,1,true,(NumericVector)P(A,_))[0];
        X = X + Dt*mu[A]*X + sigma[A]*X*norms[j]*sqrt(Dt);
        
        //check for user interrupt
        if(j % R_INTERRUPT == 0) Rcpp::checkUserInterrupt();
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
    long double mean = 0, var = 0, X;
    double CI; //NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
    {
        X = MaoSim(Dt,L,X0,p0,mu,sigma,P,E);
        
        mean = mean + (X-mean)/(j+1.0);
        if(j > 0)
            var = var + (mean-X)*(mean-X)/(j+1) - var/(j);
        
        if(j > 0 && j % N_PRINT_STATS == 0 )
        {
            CI = Rf_qnorm5(0.995,0,1,true,false)*var/sqrt(j+1);
            printf("Mao stats: %0.7Lf +- %0.5lf (DONE %.0f%%)\n",mean,CI,((float)j+1)/M*100);
        }
        
        //X[j] = MaoSim(Dt,L,X0,p0,mu,sigma,P,E);
    }
    
    return(mean);
}
