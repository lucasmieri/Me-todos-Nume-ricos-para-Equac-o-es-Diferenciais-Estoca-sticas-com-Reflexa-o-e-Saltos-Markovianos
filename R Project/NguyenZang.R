library(expm)

#MILSTEIN-TYPE PROCEDURES FOR NUMERICAL SOLUTIONS OF STOCHASTIC DIFFERENTIAL EQUATIONS WITH
#MARKOVIAN SWITCHING

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
dW = sqrt(dt)*rnorm(N)

# Euler-Maruyama:
X   = rep(0,floor(L)+1)
X[1]= X0
A   = rep(0,floor(L)+1)
A[1] = sample(c(1,2),1,prob=p0)

# TODO:
# Gerar o tempo até a primeira transição rexp(\lambd[A[1]],1)


tau <- rexp(lambda[A[1]],4)
tau<-sort(tau, decreasing = FALSE)
tau[9]=1
print(tau)
index=1
state=1


  for(j in 1:floor(L))
  {
    # TODO:
    if(tau[index] <= time[j+1]){
      index=index+1
      if(state==1){
        state=2
      }
      else {
        state=1
      }
    }
    # TROCA ESTADO
    # REGERAR O TAU como distribuição exponencial e taxa \lambda(A[j+1])
    
    #TODO: Substituir pelo Milstein
    
    
    Winc= sum(dW[(R*(j-1)+1):(R*j)])
    #X[j+1] = X[j] + Dt*mu[A[j+1]]*X[j] + sigma[A[j+1]]*X[j]*Winc
    
    X[j+1]=X[j]+Dt*X[j]*lambda[state]+X[j]*Winc*sigma[state] + 0.5*(sigma[state]^2)*X[j]*(Winc^2- Dt)
    #Xem[j] = Xem[j-1] + Dt*lamb*Xem[j-1]+ mu*Winc*Xem[j-1]
  }

# https://www.statmethods.net/advgraphs/axes.html
par(mfrow=c(1,1))
plot(time,X,main="Geometric Brownian Motion\n (Zhang's Model)",
     ylab="",xlab="Time",type="l",col="red",lwd=1)
#par(new=TRUE)
#plot(time,2-A,ylab="",xlab="Time",type="s",lty=2,col="blue",lwd=2,axes=FALSE)
#axis(4, at=c(0,1),las=1)



