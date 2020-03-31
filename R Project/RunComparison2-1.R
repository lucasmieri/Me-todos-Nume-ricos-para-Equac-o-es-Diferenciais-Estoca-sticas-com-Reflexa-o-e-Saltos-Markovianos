rm(list=ls())
library(expm)
library(Rcpp)
source("exact_mean.R")
Rcpp::sourceCpp('MaoFunction2-1.cpp')
Rcpp::sourceCpp('NguyenFunction2-1.cpp')
Rcpp::sourceCpp('TrapFunction2-1.cpp')

## Data for the simulation Zhang's Model
lambda <- c(5,10)#6.04, 8.90) #alta e baixa
sigma  <- c(0.2,0.1)#0.5,0.25) #0.44, 0.63)
prob   <- 0
r      <- c(1,-1) #1.5, -1.61)
mu     <- r + sigma^2 / 2
X0     <- 1#70.5
p0     <- c(0.5,0.5)

T  <- 1
N  <- 2^9 # number of discrete time steps
dt <- T / N
M  <- 10^7 # number of repetitions
P  <- 5     # number of step sizes

X_Mao    <- rep(0,P)
X_Nguyen <- rep(0,P)
X_Trap1  <- rep(0,P)
X_Trap2  <- rep(0,P)

load_simulations <- TRUE
if(file.exists("PlotData2-1.rda") && load_simulations)
{
    load("PlotData2-1.rda")
    
} else {
    
    for(p in 1:P)
    {
        cat("Begining for p = ",p,"...\n")
           
        R  <- 2**p
        L  <- N / R
        Dt <- R * dt
    
        # Mao's Method
        X_Mao[p]    <- MaoFuncRcpp(Dt, L, M)
    
        # Nguyen's Method
        X_Nguyen[p] <- NguyenFuncRcpp(Dt, L,M)
    
        # Trapezoidal Method
        X_Trap1[p] <- TrapFuncRcpp(1,Dt,L,M)
        X_Trap2[p] <- TrapFuncRcpp(2,Dt,L,M)
    }
    save(X_Mao,X_Nguyen,X_Trap1,X_Trap2,file="PlotData2-1.rda")
}

# Tirando a média e calculando a diferença entre médias.
RefMean   = exact_mean()
ErrMao    = abs( RefMean - X_Mao   )
ErrNguyen = abs( RefMean - X_Nguyen)
ErrTrap1  = abs( RefMean - X_Trap1 )
ErrTrap2  = abs( RefMean - X_Trap2 )

# Fazendo o gráfico
#jpeg('comparison2-1.jpg')
par(mfrow=c(2,2))
dts = c(dt*2**(1:P))

#Plot Mao
plot(dts,ErrMao,type="p",col="blue",lwd=2,log="xy",ylim = c(10^{-6},10),
     main = "Mao et al.'s Method", ylab="Log of Abs. Difference",xlab="Log of Time Step h")
reg = lm(log(ErrMao,base=10)~log(dts,base=10))
cat(sprintf("Slope Mao: %f\n",reg$coefficients[2]))
abline(reg,col="blue",lwd=3)
lines(dts,dts^(1),lwd=2,col="gray",lty=2)

#Plot Nguyen
plot(dts,ErrNguyen,type="p",col="green",lwd=2,log="xy",ylim = c(10^{-6},10),
     main = "Nguyen et al.'s Method", ylab="Log of Abs. Difference",xlab="Log of Time Step h")
reg = lm(log(ErrNguyen,base=10)~log(dts,base=10))
cat(sprintf("Slope Nguyen: %f\n",reg$coefficients[2]))
abline(reg,col="green",lwd=3)
lines(dts,dts^(1),lwd=2,col="gray",lty=2)

#Plot Trap1
plot(dts,ErrTrap1,type="p",col="red",lwd=2,log="xy",ylim = c(10^{-6},10),
     main = "Trapezoidal 1", ylab="Log of Abs. Difference",xlab="Log of Time Step h")
reg = lm(log(ErrTrap1,base=10)~log(dts,base=10))
cat(sprintf("Slope Trap 1: %f\n",reg$coefficients[2]))
abline(reg,col="red",lwd=3)
lines(dts,dts^(1),lwd=2,col="gray",lty=2)

#Plot Trap2
plot(dts,ErrTrap2,type="p",col="magenta",lwd=2,log="xy",ylim = c(10^{-6},10),
     main = "Trapezoidal 2", ylab="Log of Abs. Difference",xlab="Log of Time Step h")
reg = lm(log(ErrTrap2,base=10)~log(dts,base=10))
cat(sprintf("Slope Trap 2: %f\n",reg$coefficients[2]))
abline(reg,col="magenta",lwd=3)
lines(dts,dts^(1),lwd=2,col="gray",lty=2)

#dev.off()