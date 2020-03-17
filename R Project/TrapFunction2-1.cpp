// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

double TrapSim1(double Dt, unsigned int L, double X0, NumericVector &p0,
              NumericVector &mu, NumericVector &sigma, NumericVector &lambda, 
              IntegerVector &E)
{
    double X = X0;
    double t = 0;
    unsigned int A = sample(E,1,true,p0)[0];
    
    double theta  = 0.5;
    double alfa_1 = 1/(2*theta*(1-theta));
    double alfa_2 = ((1-theta)*(1-theta)+theta*theta)/(2*theta*(1-theta));
    
    //NumericVector norms1 = rnorm(L);
    //NumericVector norms2 = rnorm(L);
    
    double tau = rexp(1,lambda[A])[0];
    double X_step, aux;
    for(unsigned int j = 0; j < L; ++j)
    {
        t += Dt;
        if(tau <= t)
        {
            if( A == 0 ) A = 1;
            else         A = 0;
            tau += rexp(1,lambda[A])[0];
        }
        
        //step1
        X_step = X + mu[A]*X*theta*Dt + sigma[A]*X*rnorm(1)[0]*sqrt(theta*Dt);
            
        //step2
        aux = alfa_1*(X_step*sigma[A])*(X_step*sigma[A])-alfa_2*(X*sigma[A])*(X*sigma[A]);
        aux = (aux>0)? aux:0;
                
        X = X_step + (alfa_1*mu[A]*X_step - alfa_2*mu[A]*X)*(1-theta)*Dt
                + sqrt(aux)*rnorm(1)[0]*sqrt((1-theta)*Dt);

        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

double TrapSim2(double Dt, unsigned int L, double X0, NumericVector &p0,
                NumericVector &mu, NumericVector &sigma, NumericMatrix &P, IntegerVector &E)
{
    double X = X0;
    unsigned int A = sample(E,1,true,p0)[0];
    
    double theta  = 0.5;
    double alfa_1 = 1/(2*theta*(1-theta));
    double alfa_2 = ((1-theta)*(1-theta)+theta*theta)/(2*theta*(1-theta));
    
    //NumericVector norms1 = rnorm(L);
    //NumericVector norms2 = rnorm(L); 
    double X_step, aux;
    for(unsigned int j = 0; j < L; ++j)
    {
        A = sample(E,1,true,(NumericVector)P(A,_))[0];
        
        //step1
        X_step = X + mu[A]*X*theta*Dt + sigma[A]*X*rnorm(1)[0]*sqrt(theta*Dt);
            
        //step2
        aux = alfa_1*(X_step*sigma[A])*(X_step*sigma[A])-alfa_2*(X*sigma[A])*(X*sigma[A]);
        aux = (aux>0)? aux:0;
        
        X = X_step + (alfa_1*mu[A]*X_step - alfa_2*mu[A]*X)*(1-theta)*Dt
            + sqrt(aux)*rnorm(1)[0]*sqrt((1-theta)*Dt);
            
        if(j % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return(X);
}

// [[Rcpp::export]]
double TrapFuncRcpp(unsigned int type, double Dt, unsigned int L, unsigned int M)
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
    double mean = 0, alpha; //NumericVector X(M);
    if(type == 1)
        for(unsigned int j = 0; j < M; ++j)
        {
            alpha = ((double)j)/(j+1);
            mean = (alpha)*mean + (1-alpha)*TrapSim1(Dt,L,X0,p0,mu,sigma,lambda,E);
            //X[j] = TrapSim1(Dt,L,X0,p0,mu,sigma,lambda,E);
        }
    else
    {
        //Gerar P (usar expm)
        NumericMatrix Rho(2,2);
        Rho(0,0) = -lambda[0]; Rho(0,1) = lambda[0];
        Rho(1,0) = lambda[1];  Rho(1,1) = -lambda[1];
        Function expm("expm"); 
        NumericMatrix P = expm(Rho*Dt);
        for(unsigned int j = 0; j < M; ++j)
        {
            alpha = ((double)j)/(j+1);
            mean = (alpha)*mean + (1-alpha)*TrapSim2(Dt,L,X0,p0,mu,sigma,P,E);
            //X[j] = TrapSim2(Dt,L,X0,p0,mu,sigma,P,E);
        }
    }
    
    return(mean);
}
