#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector MarkovSim(double dt, unsigned int N, NumericVector p0, NumericVector lambda)
{
    IntegerVector E = {0,1};
    int A0          = sample(E,1,false,p0)[0]; 
    
    IntegerVector A(N+1);
    A[0] = A0;
    
    double t   = 0;
    double tau = rexp(1,lambda[A[0]])[0];
    
    for(unsigned int j = 0; j < N; ++j)
    {
        t += dt;
        if(tau <= t)
        {
            if( A[j] == 0 ) A[j+1] = 1;
            else            A[j+1] = 0;
            tau += rexp(1,lambda[A[j+1]])[0];
        }
        else A[j+1] = A[j];
    }
    return(A);
}

// [[Rcpp::export]]
IntegerVector CondenseMarkovSim(double R, unsigned int N, IntegerVector &v)
{
    int L = N/R;
    IntegerVector cv(L+1);
    cv[0] = v[0];
    for(unsigned int j = 1; j < L+1; ++j)
        cv[j] = v[R*j]; 
    return(cv);
}

// [[Rcpp::export]]
NumericMatrix BrownianIncs(double dt, unsigned int N, unsigned int M)
{
    NumericMatrix B(M,N);
    for(unsigned int i = 0; i < M; ++i )
        B(i, _ ) = rnorm(N)*sqrt(dt);
    return(B);
}

// [[Rcpp::export]]
NumericMatrix CondenseBrownianIncs(double R, unsigned int N, unsigned int M, NumericMatrix &v)
{
    int L = N/R;
    double sum;
    NumericMatrix cv(M,L);
    for(unsigned int k = 0; k < M; ++k)
    {
        for(unsigned int j = 0; j < L; ++j)
        {
            sum = 0.0;
            for(unsigned int i = 0; i < R; ++i)
                sum += v(k,R*j+i);
            cv(k,j) = sum;
        }
    }
    return(cv);
}

/*** R
lambda <- c(6.04, 8.90) #alta e baixa
p0     <- c(0.5,0.5)

T  <- 1.00
N  <- 2^12
dt <- T / N

t = seq(0,T,dt)
A = MarkovSim(dt,N,p0,lambda)
B = BrownianIncs(dt,N,2)

R  <- 2^5
L  <- N / R
Dt <- R * dt

tc= seq(0,T,Dt)
Ac= CondenseMarkovSim(R,N,A)
Bc= CondenseBrownianIncs(R,N,2,B)

plot(t,A,type="s")
points(tc,Ac,col="red",pch=4)

plot(tail(t,n=length(t)-1),cumsum(B[2,]),type="l")
lines(tail(tc,n=length(tc)-1),cumsum(Bc[2,]),type="l",col="red")
*/