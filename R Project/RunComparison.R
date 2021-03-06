rm(list=ls())
library(expm)
source('MaoFunction.R')
source('NguyenFunction.R')
source('TrapFunction.R')

## Data for the simulation Zhang's Model
lambda <- c(1,2)#6.04, 8.90) #alta e baixa
sigma  <- c(0.5,0.25) #0.44, 0.63)
prob   <- 0
r      <- c(1,-1) #1.5, -1.61)
mu     <- r + sigma^2 / 2
X0     <- 1#70.5
p0     <- c(0.5,0.5)

T  <- 1
N  <- 2^12 # number of discrete time steps
dt <- T / N
M  <- 1000 # number of repetitions
P  <- 5      # number of step sizes

f <- function(x) x

X_Mao    <- matrix(rep(0,M*P),ncol=P)
X_Nguyen <- matrix(rep(0,M*P),ncol=P)
X_Trap1  <- matrix(rep(0,M*P),ncol=P)
X_Trap2  <- matrix(rep(0,M*P),ncol=P)
#X_Ref    <- matrix(rep(0,M*P),ncol=1)

load_simulations <- FALSE
if(file.exists("PlotData_tmp.rda") && load_simulations)
{
    load("PlotData_tmp.rda")
    
} else {

    for( s in 1:M )
    {
        if( s%%20 == 0) 
            cat("Running form s =",s," de",M,"\n")
        
        # Run Simulations
        for(p in 1:P)
        {
            # Run Reference  (Nguyen's Method)
            #Nref <- 2^15  
            #X_Ref[(s-1)*P+p] <- tail(runNguyen(T/Nref, Nref)$X,n=1)
        
            R  <- 2**p
            L  <- N / R
            Dt <- R * dt
    
            # Mao's Method
            X_Mao[s,p] <- tail(runMao(Dt, L)$X,n=1)
    
            # Nguyen's Method
            X_Nguyen[s,p] <- tail(runNguyen(Dt, L)$X,n=1)
    
            # Trapezoidal Method
            X_Trap1[s,p] <- tail(runTrap1(Dt, L)$X,n=1)
            X_Trap2[s,p] <- tail(runTrap2(Dt, L)$X,n=1)
        }
    }

    #save(X_Ref,X_Mao,X_Nguyen,X_Trap1,X_Trap2,file="PlotData.rda")
    save(X_Mao,X_Nguyen,X_Trap1,X_Trap2,file="PlotData_tmp.rda")
}

# Tirando a média e calculando a diferença entre médias.
RefMean   = 1.6874124630113767 #mean( X_Ref[1:200000] )
ErrMao    = abs( RefMean - colMeans(X_Mao) )
ErrNguyen = abs( RefMean - colMeans(X_Nguyen) )
ErrTrap1  = abs( RefMean - colMeans(X_Trap1) )
ErrTrap2  = abs( RefMean - colMeans(X_Trap2) )

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
