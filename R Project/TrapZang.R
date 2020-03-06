library(expm)

#A WEAK TRAPEZOIDAL METHOD FOR A CLASS OF STOCHASTIC DIFFERENTIAL EQUATIONS

lambda<-c(6.04, 8.90) #alta e baixa
sigma<-c(0.44, 0.63) 
prob=0
r<-c(1.5, -1.61)
mu <- r + sigma^2/2
X0 = 70.5
p0 = c(0.5,0.5)

################## PARTE II #########################
# Fazendo os grÃ¡ficos dos caminhos

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

# TODO:
# Gerar o tempo atÃ© a primeira transiÃ§Ã£o rexp(\lambd[A[1]],1)
tau <- rexp(1,lambda[A[1]])

theta=0.5
alfa_1=1/(2*theta*(1-theta))
alfa_2=((1-theta)^2+theta^2)/(2*theta*(1-theta))


for(j in 1:floor(L))
{
    # TODO:
    if(tau <= time[j+1]){
        if( A[j] == 1 ) A[j+1] = 2
        else            A[j+1] = 1
        tau <- tau + rexp(1,lambda[A[j+1]])
    } else {
        A[j+1] = A[j]
    }

    # SUBS PELO METODO DO TRAPEZIO 
    
    #step1
    X_step=X[j]+mu[A[j+1]]*X[j]*theta*Dt+sigma[A[j+1]]*X[j]*rnorm(1)*sqrt(theta*Dt)
  
    #step2
    aux=max(0,alfa_1*(X_step*sigma[A[j+1]])^2-alfa_2*(X[j]*sigma[A[j+1]])^2) 
    
    X[j+1]=X_step+(alfa_1*mu[A[j+1]]*X_step- alfa_2*mu[A[j+1]]*X[j])*(1-theta)*Dt
           + sqrt(aux)*rnorm(1)*sqrt((1-theta)*Dt)
}

# https://www.statmethods.net/advgraphs/axes.html
par(mfrow=c(1,1))
plot(time,X,main="Geometric Brownian Motion\n (Zhang's Model)",
     ylab="",xlab="Time",type="l",col="red",lwd=1)
par(new=TRUE)
plot(time,2-A,ylab="",xlab="Time",type="s",lty=2,col="blue",lwd=2,axes=FALSE)
axis(4, at=c(0,1),las=1)



