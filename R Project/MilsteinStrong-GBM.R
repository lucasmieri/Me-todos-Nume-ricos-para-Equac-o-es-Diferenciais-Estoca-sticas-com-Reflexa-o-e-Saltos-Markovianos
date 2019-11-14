lamb  = 2
mu    = 1
Xzero = 1

################## PARTE I #########################
# Fazendo os gráficos dos erros
T  = 1
N  = 2**9 # number of discrete time steps
dt = T/N
M  = 2000 # number of repetitions
P  = 5
Xerr   = matrix(rep(0,M*P),ncol=P)

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
            Xem = Xem + lamb*Xem*Dt + mu*Xem*Winc + 0.5*(mu^2)*Xem*(Winc^2-Dt)
        }
        Xerr[s,p]=abs(Xem-Xtrue)
    }
}

# Tirando a média do erro para cada p 
mXerr = colMeans(Xerr)

# Fazendo o gráfico do erro
dts = c(dt*2**(1:5))
plot(dts,mXerr,ylab="Log of Abs. Difference",
     xlab="Log of Time Step h",type="b",col="blue",ylim=c(0.001,10),lwd=2,log="xy")
lines(dts,dts,lwd=2,col="red")
legend(0.005, 5, legend=c("Milstein", "Reference (slope=1)"),
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
    W = sum(dW[0:(R*j)])
    XtrueV[j+1] = Xzero*exp((lamb-0.5*mu**2)*(j*Dt)+mu*W)
}
    
# Euler-Maruyama:
XemV = rep(0,floor(L)+1)
XemV[1]=Xzero
for(j in 1:floor(L))
{
    Winc= sum(dW[(R*(j-1)+1):(R*j)])
    XemV[j+1] = XemV[j] + lamb*XemV[j]*Dt + mu*XemV[j]*Winc + 0.5*(mu^2)*XemV[j]*(Winc^2-Dt)
}
plot(time,XtrueV,ylab="",xlab="Time",type="l",col="red",ylim=c(0,10))
lines(time,XemV,type="l",lwd=2,col="blue")
legend(0.005,9, legend=c("Euler-Maruyama", "Analytic Solution"),
       col=c("blue", "red"), lty=1, lwd=2)
