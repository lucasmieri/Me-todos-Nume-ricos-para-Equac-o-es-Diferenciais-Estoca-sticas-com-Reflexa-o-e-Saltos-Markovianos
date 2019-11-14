theta = 2
mu    = 0.5
Xzero = 1
sigma = 5

################## PARTE I #########################
# Fazendo os gráficos dos erros
T  = 1
N  = 2**9 # number of discrete time steps
dt = T/N
M  = 2000 # number of repetitions
P  = 5
Xerr = matrix(rep(0,M*P),ncol=P)
for(s in 1:M)
{
    dW = sqrt(dt)*rnorm(N)
    WT = sum(sigma*exp(-theta*(T-seq(dt,T,dt))) %*% dW)
    
    # Analytic solution
    Xtrue = Xzero*exp(-theta*T) + mu*(1-exp(-theta*T)) + WT

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
        Xerr[s,p]=abs(Xem-Xtrue)
    }
}

# Tirando a média do erro para cada p 
mXerr = colMeans(Xerr)

# Fazendo o gráfico
dts = c(dt*2**(1:P))
plot(dts,mXerr,main="Strong Convergence Ornstein Uhlenbeck",
     ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",lwd=2,ylim=c(0.001,10),log="xy")
lines(dts,dts^1,lwd=2,col="red")
legend(0.005, 1e+01, legend=c("Euler-Maruyama", "Reference (slope=1)"),
       col=c("blue", "red"), lty=1, lwd=2)

################## PARTE II #########################
# Fazendo os gráficos dos caminhos
p = 1
R=2**p
L=N/R
Dt=R*dt

set.seed(15) # Fazendo a semente ser 15 para dar sempre o mesmo caminho
time = seq(0,T,Dt)
dW = sqrt(dt)*rnorm(N)

# Analytic solution
XtrueV = rep(0,floor(L)+1)
XtrueV[1] = Xzero
for(j in 1:floor(L)){
    t = j*Dt
    WT = sum(sigma*exp(-theta*(t-seq(dt,t,dt))) %*% dW[0:(R*j)])
    XtrueV[j+1] = Xzero*exp(-theta*t) + mu*(1-exp(-theta*t)) + WT
}

# Euler-Maruyama:
XemV = rep(0,floor(L)+1)
XemV[1]=Xzero
for(j in 1:floor(L))
{
    Winc= sum(dW[(R*(j-1)+1):(R*j)])
    XemV[j+1] = XemV[j] + Dt*theta*(mu-XemV[j]) + sigma*Winc
}
plot(time,XtrueV,main="Ornstein Uhlenbeck",
    ylab="",xlab="Time",type="l",col="red",ylim=c(-5,10))
lines(time,XemV,type="l",lwd=2,col="blue")
legend(0.005,10, legend=c("Euler-Maruyama", "Analytic Solution"),
       col=c("blue", "red"), lty=1, lwd=2)
