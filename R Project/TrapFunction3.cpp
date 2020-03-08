#include <Rcpp.h>
using namespace Rcpp;

double TrapSim(double Dt, unsigned int L, double X0, 
              NumericVector mu, NumericVector sigma,
              IntegerVector A, NumericVector B)
{
    double X = X0;
    
    double theta  = 0.5;
    double alfa_1 = 1/(2*theta*(1-theta));
    double alfa_2 = ((1-theta)*(1-theta)+theta*theta)/(2*theta*(1-theta));
   
    NumericVector norms = rnorm(L)*sqrt(Dt); 
    double X_step, aux;
    for(unsigned int j = 1; j < L+1; ++j)
    {
        //step1
        X_step = X + mu[A[j]]*X*theta*Dt + sigma[A[j]]*X*norms[j-1]*sqrt(theta);
            
        //step2
        aux = alfa_1*(X_step*sigma[A[j]])*(X_step*sigma[A[j]])-alfa_2*(X*sigma[A[j]])*(X*sigma[A[j]]);
        aux = (aux>0)? aux:0;
                
        X = X_step + (alfa_1*mu[A[j]]*X_step - alfa_2*mu[A[j]]*X)*(1-theta)*Dt
                + sqrt(aux)*B[j-1]*sqrt((1-theta));

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
NumericVector TrapFuncRcpp(double Dt, unsigned int L, unsigned int M,
                           IntegerVector A, NumericMatrix B)
{
    //Global env data
    Environment env = Environment::global_env();
    double X0            = env["X0"];
    NumericVector mu     = env["mu"];
    NumericVector sigma  = env["sigma"];
    
    //Running Simulations 
    NumericVector X(M);
    for(unsigned int j = 0; j < M; ++j)
        X[j] = TrapSim(Dt,L,X0,mu,sigma,A,B(j,_));
    
    return(X);
}
