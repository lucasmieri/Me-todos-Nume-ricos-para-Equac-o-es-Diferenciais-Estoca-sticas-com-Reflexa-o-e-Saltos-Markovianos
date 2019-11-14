theta = 2
mu    = 0.5
Xzero = 1
sigma = 5

################## PARTE III #########################
# Fazendo os gráficos do erro das médias (FRACO)
T  = 1
N  = 2**9 # number of discrete time steps
dt = T/N
M  = 2000 # number of repetitions
P  = 5

Xem_at_T   = matrix(rep(0,M*P),ncol=P)
Xtrue_at_T = rep(0,M)

for(s in 1:M)
{
    dW = sqrt(dt)*rnorm(N)
    WT = sum(sigma*exp(-theta*(T-seq(dt,T,dt))) %*% dW)
    
    # Analytic solution
    Xtrue = Xzero*exp(-theta*T) + mu*(1-exp(-theta*T)) + WT
    Xtrue_at_T[s] = Xtrue

    # Euler-Maruyama:
    for( p in 1:P )
    {
        R=2**p
        L=N/R
        Dt=R*dt
        Xem=Xzero
        for( j in 1:floor(L))
        {
            Winc= sum(dW[(R*(j-1)+1):(R*j)])
            Xem = Xem + theta*(mu-Xem)*Dt + sigma*Winc
        }
        Xem_at_T[s,p]=Xem
    }
}

# Tirando a média e calculando a diferença entre médias.
mXerr = abs(mean(Xtrue_at_T) - colMeans(Xem_at_T))

# Fazendo o gráfico
dts = c(dt*2**(1:P))
plot(dts,mXerr,main="Weak Convergence Ornstein Uhlenbeck", ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,ylim=c(0.0001,10),log="xy")
lines(dts,dts^1,lwd=2,col="red")
legend(0.005, 1e+01, legend=c("Euler-Maruyama", "Reference (slope=1)"),
       col=c("blue", "red"), lty=1, lwd=2)
