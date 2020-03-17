// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
#define R_INTERRUPT 1000
#define N_PRINT_STATS 1000

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
    
    NumericVector norms1 = rnorm(L);
    NumericVector norms2 = rnorm(L);
    
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
        X_step = X + mu[A]*X*theta*Dt + sigma[A]*X*norms2[j]*sqrt(theta*Dt);
            
        //step2
        aux = alfa_1*(X_step*sigma[A])*(X_step*sigma[A])-alfa_2*(X*sigma[A])*(X*sigma[A]);
        aux = (aux>0)? aux:0;
                
        X = X_step + (alfa_1*mu[A]*X_step - alfa_2*mu[A]*X)*(1-theta)*Dt
                + sqrt(aux)*norms1[j]*sqrt((1-theta)*Dt);

        if(j % R_INTERRUPT == 0) Rcpp::checkUserInterrupt();
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
    
    NumericVector norms1 = rnorm(L);
    NumericVector norms2 = rnorm(L); 
    NumericVector unifs  = runif(L);
    double X_step, aux;
    for(unsigned int j = 0; j < L; ++j)
    {
        if(unifs[j] <= P(A,0)) A = 0;
        else                   A = 1;
        //A = sample(E,1,true,(NumericVector)P(A,_))[0];
        
        //step1
        X_step = X + mu[A]*X*theta*Dt + sigma[A]*X*norms2[j]*sqrt(theta*Dt);
            
        //step2
        aux = alfa_1*(X_step*sigma[A])*(X_step*sigma[A])-alfa_2*(X*sigma[A])*(X*sigma[A]);
        aux = (aux>0)? aux:0;
        
        X = X_step + (alfa_1*mu[A]*X_step - alfa_2*mu[A]*X)*(1-theta)*Dt
            + sqrt(aux)*norms1[j]*sqrt((1-theta)*Dt);
            
        if(j % R_INTERRUPT == 0)  Rcpp::checkUserInterrupt();
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
    double mean = 0, var = 0, X, CI; //NumericVector X(M);
    if(type == 1)
        for(unsigned int j = 0; j < M; ++j)
        {
            //X[j] = TrapSim1(Dt,L,X0,p0,mu,sigma,lambda,E);
            X = TrapSim1(Dt,L,X0,p0,mu,sigma,lambda,E);
            
            mean = mean + (X-mean)/(j+1);
            if(j > 0)
                var = var + (mean-X)*(mean-X)/(j+1) - var/(j);
            
            if(j > 0 && j % N_PRINT_STATS == 0 )
            {
                CI = Rf_qnorm5(0.995,0,1,true,false)*var/sqrt(j+1);
                printf("Trap1 stats: %0.7f +- %0.5f (DONE %.2f%%)\n",mean,CI,((float)j+1)/M);
            }
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
            //X[j] = TrapSim2(Dt,L,X0,p0,mu,sigma,P,E);
            X = TrapSim2(Dt,L,X0,p0,mu,sigma,P,E);
            
            mean = mean + (X-mean)/(j+1);
            if(j > 0)
                var = var + (mean-X)*(mean-X)/(j+1) - var/(j);
            
            if(j > 0 && j % N_PRINT_STATS == 0 )
            {
                CI = Rf_qnorm5(0.995,0,1,true,false)*var/sqrt(j+1);
                printf("Trap2 stats: %0.7f +- %0.5f (DONE %.2f%%)\n",mean,CI,((float)j+1)/M);
            }
        }
    }
    
    return(mean);
}
