lamb  = 2
mu    = 1
Xzero = 1

### Definir os parâmetos usando os dados do artigo do Zang
# TODO
#lambda_1 = 
#lambda_2 = 
#mu_1 = 
#mu_2 = 
#sigma_1 = 
#sigma_2 = 

################## PARTE II #########################
# Fazendo os gráficos dos caminhos
p = 1
R=2**p
L=N/R
Dt=R*dt

set.seed(15) # Fazendo a semente ser 15 para dar sempre o mesmo caminho
time = seq(0,T,Dt)
dW = sqrt(dt)*rnorm(N)

# Euler-Maruyama:
XemV = rep(0,floor(L)+1)
XemV[1]=Xzero
for(j in 1:floor(L))
{
    # TODO
    #Gerar dinamica da cadeia de Markov
    # 1) Gerar P (usar expm)
    # 2) Gerar o estado, gerando um número aleatório uniforme em (0,1)
    
    Winc= sum(dW[(R*(j-1)+1):(R*j)])
    
    #Usar os paramatros dependente da cadeia de Markov
    # TODO
    XemV[j+1] = XemV[j] + Dt*lamb*XemV[j] + mu*XemV[j]*Winc
}
plot(time,XtrueV,main="Geometric Brownian Motion\n (Zhang's Model)",
     ylab="",xlab="Time",type="l",col="red",ylim=c(0,10))