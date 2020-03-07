lamb  = 2
mu    = 1
Xzero = 1

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
    WT  = sum(dW)
    
    # Analytic solution
    Xtrue = Xzero*exp((lamb-0.5*(mu)**2)*T+mu*WT)
    Xtrue_at_T[s] = Xtrue
    
    # Euler-Maruyama:
    for( p in 1:5 )
    {
        R=2**p
        L=N/R
        Dt=R*dt
        Xem=Xzero
        for( j in 1:floor(L))
        {
            Winc= sum(dW[(R*(j-1)+1):(R*j)])
            Xem = Xem + lamb*Xem*Dt + mu*Xem*Winc + 0.5*(mu^2)*Xem*(Winc^2-Dt)
        }
        Xem_at_T[s,p]=Xem
    }
}

# Tirando a média e calculando a diferença entre médias.
mXerr = abs(mean(Xtrue_at_T) - colMeans(Xem_at_T))

# Fazendo o gráfico do erro
dts = c(dt*2**(1:5))
plot(dts,mXerr,main="Weak Convergence Geometric Brownian Motion",ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",ylim=c(0.001,10),lwd=2,log="xy")
lines(dts,dts,lwd=2,col="red")
legend(0.005, 5, legend=c("Milstein", "Reference (slope=1)"),
       col=c("blue", "red"), lty=1, lwd=2)