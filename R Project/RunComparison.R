rm(list=ls())
library(expm)
source('MaoZangFunction.R')

## Data for the simulation
lambda <- c(6.04, 8.90) #alta e baixa
sigma  <- c(0.44, 0.63) 
prob   <- 0
r      <- c(1.5, -1.61)
mu     <- r + sigma^2/2
X0     <- 70.5
p0     <- c(0.5,0.5)

T  <- 1
N  <- 2**12 # number of discrete time steps
dt <- T/N
M  <- 2000 # number of repetitions
P  <- 5    # number of step sizes

f <- function(x) x

Fval_paths_MaoZang <- matrix(rep(0,M*P),ncol=P)
for( s in 1:M )
{
    for( p in 1:P )
    {
        R  <- 2**p
        L  <- N/R
        Dt <- R*dt
        
        X <- runMaoZang(Dt,L)$X
        Fval_paths_MaoZang[s,p] <- f(X[length(X)])
    }
}