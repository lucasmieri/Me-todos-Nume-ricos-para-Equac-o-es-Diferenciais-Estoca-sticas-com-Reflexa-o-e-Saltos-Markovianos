runTrap1 <- function(Dt,L)
{
    theta=0.5
    alfa_1=1/(2*theta*(1-theta))
    alfa_2=((1-theta)^2+theta^2)/(2*theta*(1-theta))
    
    # Initializing
    X    = rep(0,floor(L)+1)
    X[1] = X0
    A    = rep(0,floor(L)+1)
    A[1] = sample(c(1,2),1,prob=p0)

    time = seq(0,T,Dt)
    tau  = rexp(1,lambda[A[1]])
    
    for(j in 1:floor(L))
    {
        # TODO:
        if(tau <= time[j+1]){
            if( A[j] == 1 ) A[j+1] = 2
            else            A[j+1] = 1
            tau = tau + rexp(1,lambda[A[j+1]])
        } else {
            A[j+1] = A[j]
        }
        
        #step1
        X_step = X[j] + mu[A[j+1]]*X[j]*theta*Dt + sigma[A[j+1]]*X[j]*rnorm(1)*sqrt(theta*Dt)
        
        #step2
        aux = max(0,alfa_1*(X_step*sigma[A[j+1]])^2-alfa_2*(X[j]*sigma[A[j+1]])^2)
        
        X[j+1] = X_step + (alfa_1*mu[A[j+1]]*X_step - alfa_2*mu[A[j+1]]*X[j])*(1-theta)*Dt
        + sqrt(aux)*rnorm(1)*sqrt((1-theta)*Dt)
    }
    return(list('X' = X, 'A' = A))
}

runTrap2 <- function(Dt,L)
{
    theta=0.5
    alfa_1=1/(2*theta*(1-theta))
    alfa_2=((1-theta)^2+theta^2)/(2*theta*(1-theta))
    
    # Initializing
    X    = rep(0,floor(L)+1)
    X[1] = X0
    A    = rep(0,floor(L)+1)
    A[1] = sample(c(1,2),1,prob=p0)
    
    # Gerar P (usar expm)
    Rho <- matrix(c(-lambda[1], lambda[1], lambda[2], -lambda[2]), nrow=2,byrow=TRUE)
    P<-expm(Rho*Dt)
    
    for(j in 1:floor(L))
    {
        # Gerar dinamica da cadeia de Markov
        A[j+1] = sample(c(1,2),1,prob=P[A[j],]) 

        #step1
        X_step = X[j] + mu[A[j+1]]*X[j]*theta*Dt + sigma[A[j+1]]*X[j]*rnorm(1)*sqrt(theta*Dt)
        
        #step2
        aux = max(0,alfa_1*(X_step*sigma[A[j+1]])^2-alfa_2*(X[j]*sigma[A[j+1]])^2) 
        
        X[j+1] = X_step + (alfa_1*mu[A[j+1]]*X_step- alfa_2*mu[A[j+1]]*X[j])*(1-theta)*Dt
        + sqrt(aux)*rnorm(1)*sqrt((1-theta)*Dt)
    }
    return(list('X' = X, 'A' = A))
}

