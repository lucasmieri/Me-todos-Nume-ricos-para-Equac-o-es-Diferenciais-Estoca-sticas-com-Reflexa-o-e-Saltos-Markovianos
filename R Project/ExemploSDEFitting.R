library("Sim.DiffProc")

# Gerando dados da SDE a partir de uma simulação
f <- expression( 1*x )
g <- expression( 2*x )
sim <- snssde1d(drift=f,diffusion=g,x0=2,N=10^4,Dt=10^-4)
mydata <- sim$X

# Fazendo o fitting usando a biblioteca Sim.DiffProc
fx <- expression( theta[1]*x ) ## drift coefficient of model
gx <- expression( theta[2]*x ) ## diffusion coefficient of model 
fitmod <- fitsde(data = mydata, drift = fx, 
           diffusion = gx, 
           start = list(theta1=1, theta2=1),
           pmle="euler")

fitmod
print(coef(fitmod))