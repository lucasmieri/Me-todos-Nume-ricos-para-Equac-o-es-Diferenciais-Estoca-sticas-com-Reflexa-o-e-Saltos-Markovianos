rm(list=ls())
library(expm)
library(Rcpp)
Rcpp::sourceCpp('MaoFunction3.cpp')
Rcpp::sourceCpp('NguyenFunction3.cpp')
Rcpp::sourceCpp('TrapFunction3.cpp')
Rcpp::sourceCpp('Utils.cpp')

## Data for the simulation Zhang's Model
lambda <- c(6.04, 8.90) #alta e baixa
sigma  <- c(0.44, 0.63)
prob   <- 0
r      <- c(1.5, -1.61)
mu     <- r + sigma^2 / 2
X0     <- 70.5
p0     <- c(0.5,0.5)

T  <- 1.00
N  <- 2^17 # number of discrete time steps
dt <- T / N
M  <- 2000 # number of repetitions
C  <- 200  # number of markov chains
P  <- 5    # number of step sizes

f <- function(x) x

X_Mao    <- matrix(rep(0,M*P*C),ncol=P)
X_Nguyen <- matrix(rep(0,M*P*C),ncol=P)
X_Trap   <- matrix(rep(0,M*P*C),ncol=P)
X_Ref    <- matrix(rep(0,M*C),ncol=1)

run_simulations <- TRUE
if(file.exists("PlotData3.rda"))
{
    load("PlotData3.rda")
    run_simulations <- FALSE
    
} else {

    # Run Simulations
    for(c in 1:C)
    {
        cat("Begining for c =",c,"\n")
        # Gen Markov chain
        A = MarkovSim(dt,N,p0,lambda)
        # Gen Browninan Motion
        B = BrownianIncs(dt,N,M)
        
        # Run Reference  (Nguyen's Method)
        X_Ref[((c-1)*M+1):(c*M)] <- NguyenFuncRcpp(dt,N,M,A,B)
        
        cat("  done with reference...\n") 
    
        for(p in 1:P)
        {
            cat("  p =",p,"\n") 
            
            R  <- 2**(p+2)
            L  <- N / R
            Dt <- R * dt
            
            # Condense
            Ac = CondenseMarkovSim(R,N,A)
            rm(A)
            Bc = CondenseBrownianIncs(R,N,M,B)
            rm(B)
        
            # Mao's Method
            X_Mao[((c-1)*M+1):(c*M),p]    <- MaoFuncRcpp(Dt,L,M,Ac,Bc)
        
            # Nguyen's Method
            X_Nguyen[((c-1)*M+1):(c*M),p] <- NguyenFuncRcpp(Dt,L,M,Ac,Bc)
        
            # Trapezoidal Method
            X_Trap[((c-1)*M+1):(c*M),p]   <- TrapFuncRcpp(Dt,L,M,Ac,Bc)
        }
    }
    save(X_Ref,X_Mao,X_Nguyen,X_Trap,file="PlotData3.rda")
}

# Tirando a média e calculando a diferença entre médias.
RefMean   = mean( X_Ref )   #apply(X_Ref,1:2,f) )
ErrMao    = abs( RefMean - colMeans(X_Mao)    ) #apply(X_Mao,   1:2,f)) )
ErrNguyen = abs( RefMean - colMeans(X_Nguyen) ) #apply(X_Nguyen,1:2,f)) )
ErrTrap   = abs( RefMean - colMeans(X_Trap)   ) #apply(X_Trap1, 1:2,f)) )

# Fazendo o gráfico
dts = c(dt*2**(1:P))
plot(dts,ErrMao,main="Weak Convergence For Zhang's Model",
     ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,log="xy",ylim=c(0.0001,20))
lines(dts,ErrTrap,lwd=2,col="red",type="b")
lines(dts,ErrNguyen,lwd=2,col="green",type="b")
lines(dts,10*dts^(1),lwd=2,col="cyan",lty=2)
#legend(0.0022, 25,
#       legend=c("Mao",
#                "Nguyen","Trap1", "Trap2", "Reference (slope=1)"),
#       col=c("blue", "red", "magenta", "green","cyan"), lty=1, lwd=2)
