lamb  = 2
mu    = 1
Xzero = 1

T  = 1
N  = 2**9 # number of discrete time steps
dt = T/N
M  = 1000 # number of repetitions
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
            Xem = Xem + Dt*lamb*Xem + mu*Xem*Winc
        }
        Xerr[s,p]=abs(Xem-Xtrue)
    }
}

# Tirando a média do erro para cada p 
mXerr = colMeans(Xerr)

# Fazendo o gráfico
dts = c(dt*2**(1:5))
plot(dts,mXerr,type="l",col="blue",ylim=c(0.01,10),log="xy")
lines(dts,dts^(1/2),col="red")
