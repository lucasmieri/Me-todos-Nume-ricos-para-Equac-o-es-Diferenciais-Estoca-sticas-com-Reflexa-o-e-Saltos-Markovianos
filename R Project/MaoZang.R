library(expm)

lambda<-c(6.04, 8.90) #alta e baixa
sigma<-c(0.44, 0.63) 
prob=0
r<-c(1.5, -1.61)
mu <- r + sigma^2/2
X0 = 70.5
p0 = c(0.5,0.5)

################## PARTE II #########################
# Fazendo os gráficos dos caminhos

N= 2**12 # number of discrete time steps
R=2
L=N/R
T=1
dt = T/N
Dt=R*dt

#set.seed(15) # Fazendo a semente ser 15 para dar sempre o mesmo caminho
time = seq(0,T,Dt)

# Euler-Maruyama:
X   = rep(0,floor(L)+1)
X[1]= X0
A   = rep(0,floor(L)+1)
A[1] = sample(c(1,2),1,prob=p0)

# Gerar P (usar expm)
Rho <- matrix(c(-lambda[1], lambda[1], lambda[2], -lambda[2]), nrow=2,byrow=TRUE)
P<-expm(Rho*Dt)

for(j in 1:floor(L))
{
    # Gerar dinamica da cadeia de Markov
    A[j+1] = sample(c(1,2),1,prob=P[A[j],])

    Winc= sqrt(Dt)*rnorm(1)
    X[j+1] = X[j] + Dt*mu[A[j+1]]*X[j] + sigma[A[j+1]]*X[j]*Winc
}

# https://www.statmethods.net/advgraphs/axes.html
par(mfrow=c(1,1))
plot(time,X,main="Geometric Brownian Motion\n (Zhang's Model)",
     ylab="",xlab="Time",type="l",col="red",lwd=1)
par(new=TRUE)
plot(time,2-A,ylab="",xlab="Time",type="s",lty=2,col="blue",lwd=2,axes=FALSE)
axis(4, at=c(0,1),las=1)

