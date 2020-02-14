### Definir os parâmetos usando os dados do artigo do Zang
# TODO
#lambda_1 = 
#lambda_2 = 
#mu_1 = 
#mu_2 = 
#sigma_1 = 
#sigma_2 = 

lambda<-c(6.04, 8.90) #alta e baixa
sigma<-c(0.44, 0.63) 
prob=0
r<-c(1.5, -1.61)
p<- matrix(c(-lambda[1], lambda[1], lambda[2], -lambda[2]), nrow=2)
Q<-expm(p)
Xzero = 70.5

################## PARTE II #########################
# Fazendo os gráficos dos caminhos

N= 2**9 # number of discrete time steps
R=2
L=N/R
T=1
dt = T/N
Dt=R*dt



set.seed(15) # Fazendo a semente ser 15 para dar sempre o mesmo caminho
time = seq(0,T,Dt)
dW = sqrt(dt)*rnorm(N)

# Euler-Maruyama:
XemV = rep(0,floor(L)+1)
XemV[1]=Xzero
Xemv_zhang<-XemV


#optimal solution
#dx= rdt+ sigma dw

for(j in 1:floor(L))
{
    
    # TODO
    #Gerar dinamica da cadeia de Markov
    # 1) Gerar P (usar expm)
    
    
    # 2) Gerar o estado, gerando um número aleatório uniforme em (0,1)
    
    Winc= sum(dW[(R*(j-1)+1):(R*j)])
    
    #Usar os paramatros dependente da cadeia de Markov
    # TODO
    
    
    
    
   # XemV[j+1] = XemV[j] + Dt*lambda[prob+1]*XemV[j] + sigma[prob+1]*XemV[j]*Winc
    
    #zhang
    
    Xemv_zhang[j+1] = Xemv_zhang[j] + Dt*r[prob+1]*XemV[j] + sigma[prob+1]*XemV[j]*Winc
    
    prob<- sample(0:1, 1)
    
}




plot(time,Xemv_zhang,main="Geometric Brownian Motion\n (Zhang's Model)",
     ylab="",xlab="Time",type="l",col="red")
