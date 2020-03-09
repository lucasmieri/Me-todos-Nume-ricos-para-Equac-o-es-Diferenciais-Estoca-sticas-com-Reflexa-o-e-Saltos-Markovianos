rm(list=ls())
library(expm)
library(Rcpp)
Rcpp::sourceCpp('MaoFunction2.cpp')
Rcpp::sourceCpp('NguyenFunction2.cpp')
Rcpp::sourceCpp('TrapFunction2.cpp')

## Data for the simulation Zhang's Model
lambda <- c(6.04, 8.90) #alta e baixa
sigma  <- c(0.44, 0.63)
prob   <- 0
r      <- c(1.5, -1.61)
mu     <- r + sigma^2 / 2
X0     <- 70.5
p0     <- c(0.5,0.5)

T  <- 1.00
N  <- 2^12 # number of discrete time steps
dt <- T / N
M  <- 10000000 # number of repetitions
P  <- 5        # number of step sizes

f <- function(x) x

X_Mao    <- matrix(rep(0,M*P),ncol=P)
X_Nguyen <- matrix(rep(0,M*P),ncol=P)
X_Trap1  <- matrix(rep(0,M*P),ncol=P)
X_Trap2  <- matrix(rep(0,M*P),ncol=P)
X_Ref    <- matrix(rep(0,M*1),ncol=1)

run_simulations <- TRUE
if(file.exists("PlotData2.rda"))
{
    load("PlotData2.rda")
    run_simulations <- FALSE
    
} else {

    # Run Reference  (Nguyen's Method)
    Nref  <- 2^15
    X_Ref <- NguyenFuncRcpp(T/Nref, Nref, M)
    cat("Done with reference...\n") 
    
    for(p in 1:P)
    {
        cat("Begining for p = ",p,"...\n")
           
        R  <- 2**p
        L  <- N / R
        Dt <- R * dt
    
        # Mao's Method
        X_Mao[,p]    <- MaoFuncRcpp(Dt, L, M)
    
        # Nguyen's Method
        X_Nguyen[,p] <- NguyenFuncRcpp(Dt, L,M)
    
        # Trapezoidal Method
        X_Trap1[,p] <- TrapFuncRcpp(1,Dt,L,M)
        X_Trap2[,p] <- TrapFuncRcpp(2,Dt,L,M)
    }
    save(X_Ref,X_Mao,X_Nguyen,X_Trap1,X_Trap2,file="PlotData2.rda")
}

# Tirando a média e calculando a diferença entre médias.
RefMean   = mean( X_Ref )   #apply(X_Ref,1:2,f) )
ErrMao    = abs( RefMean - colMeans(X_Mao)    ) #apply(X_Mao,   1:2,f)) )
ErrNguyen = abs( RefMean - colMeans(X_Nguyen) ) #apply(X_Nguyen,1:2,f)) )
ErrTrap1  = abs( RefMean - colMeans(X_Trap1)  ) #apply(X_Trap1, 1:2,f)) )
ErrTrap2  = abs( RefMean - colMeans(X_Trap2)  ) #apply(X_Trap2, 1:2,f)) )

# Fazendo o gráfico
jpeg('comparison2.jpg')
dts = c(dt*2**(1:P))
plot(dts,ErrMao,main="Weak Convergence For Zhang's Model",
     ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,ylim=c(0.0001,20),log="xy")
lines(dts,ErrTrap1,lwd=2,col="red",type="b")
lines(dts,ErrTrap2,lwd=2,col="magenta",type="b")
lines(dts,ErrNguyen,lwd=2,col="green",type="b")
lines(dts,dts^(1),lwd=2,col="cyan",lty=2)
dev.off()
#legend(0.0022, 25,
#       legend=c("Mao",
#                "Nguyen","Trap1", "Trap2", "Reference (slope=1)"),
#       col=c("blue", "red", "magenta", "green","cyan"), lty=1, lwd=2)
