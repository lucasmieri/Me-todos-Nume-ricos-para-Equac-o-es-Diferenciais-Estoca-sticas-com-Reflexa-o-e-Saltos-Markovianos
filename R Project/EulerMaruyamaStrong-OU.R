theta = 2
mu    = 1
Xzero = 1
sigma = 2

T  = 1
N  = 2**9 # number of discrete time steps
dt = T/N
M  = 2000 # number of repetitions
P  = 5
Xerr = matrix(rep(0,M*P),ncol=P)

for(s in 1:M)
{
    dW = sqrt(dt)*rnorm(N)
    W  = cumsum(dW)
    
    # Analytic solution
    Xtrue = Xzero*exp((lamb-0.5*mu**2)*T+mu*W[length(W)])
    
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
            Xem = Xem + Dt*theta*(mu-Xem) + sigma*Winc
        }
        Xerr[s,p]=abs(Xem-Xtrue)
    }
}

# Tirando a média do erro para cada p 
mXerr = colMeans(Xerr)

# Fazendo o gráfico
dts = c(dt*2**(1:5))
plot(dts,mXerr,ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,ylim=c(0.001,10),log="xy")
lines(dts,dts^(1/2),lwd=2,col="red")
legend(0.005, 5, legend=c("Euler-Maruyama", "Reference (slope=1/2)"),
       col=c("blue", "red"), lty=1, lwd=2)
