rm(list=ls())
library(expm)
source('MaoFunction.R')
source('NguyenFunction.R')
source('TrapFunction.R')

## Data for the simulation Zhang's Model
lambda <- c(6.04, 8.90) #alta e baixa
sigma  <- c(0.44, 0.63)
prob   <- 0
r      <- c(1.5, -1.61)
mu     <- r + sigma^2 / 2
X0     <- 70.5
p0     <- c(0.5,0.5)

T  <- 1
N  <- 2^12 # number of discrete time steps
dt <- T / N
M  <- 2000 # number of repetitions
P  <- 5    # number of step sizes

f <- function(x) x

Fval_paths_Mao    <- matrix(rep(0,M*P),ncol=P)
Fval_paths_Nguyen <- matrix(rep(0,M*P),ncol=P)
Fval_paths_Trap1  <- matrix(rep(0,M*P),ncol=P)
Fval_paths_Trap2  <- matrix(rep(0,M*P),ncol=P)
Fval_paths_REF    <- matrix(rep(0,M),ncol=1)

if(file.exists("PlotData.rda"))
{
    load("PlotData.rda")
    
} else {
    
    for( s in 1:M )
    {
        # Run Reference
        # Nguyen's Method
        Nref <- 2^15  
        X <- runNguyen(T/Nref, Nref)$X
        Fval_paths_REF[s,1] <- f(X[length(X)])
        
        # Run Simulations
        for(p in 1:P)
        {
            R  <- 2**p
            L  <- N / R
            Dt <- R * dt
    
            # Mao's Method
            X <- runMao(Dt, L)$X
            Fval_paths_Mao[s, p] <- f(X[length(X)])
    
            # Nguyen's Method
            X <- runNguyen(Dt, L)$X
            Fval_paths_Nguyen[s, p] <- f(X[length(X)])
    
            # Trapezoidal Method
            X <- runTrap1(Dt, L)$X
            Fval_paths_Trap1[s, p] <- f(X[length(X)])
    
            X <- runTrap2(Dt, L)$X
            Fval_paths_Trap2[s, p] <- f(X[length(X)])
        }
    }
    save(Fval_paths_REF,
         Fval_paths_Mao,
         Fval_paths_Nguyen,
         Fval_paths_Trap1,
         Fval_paths_Trap2,file="PlotData.rda")
}

# Tirando a média e calculando a diferença entre médias.
RefMean   = mean(Fval_paths_REF)
ErrMao    = abs(RefMean - colMeans(Fval_paths_Mao))
ErrNguyen = abs(RefMean - colMeans(Fval_paths_Nguyen))
ErrTrap1  = abs(RefMean - colMeans(Fval_paths_Trap1))
ErrTrap2  = abs(RefMean - colMeans(Fval_paths_Trap2))

# Fazendo o gráfico
dts = c(dt*2**(1:5))
plot(dts,ErrMao,main="Weak Convergence For Zhang's Model",
     ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,ylim=c(0.0001,20),log="xy")
lines(dts,ErrTrap1,lwd=2,col="red",type="b")
lines(dts,ErrTrap2,lwd=2,col="magenta",type="b")
lines(dts,ErrNguyen,lwd=2,col="green",type="b")
lines(dts,dts^(1),lwd=2,col="cyan",lty=2)
#legend(0.0022, 25,
#       legend=c("Mao",
#                "Nguyen","Trap1", "Trap2", "Reference (slope=1)"),
#       col=c("blue", "red", "magenta", "green","cyan"), lty=1, lwd=2)