runMao <- function(Dt,L)
{
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
        A[j+1] <- sample(c(1,2),1,prob=P[A[j],])

        Winc   <- sqrt(Dt)*rnorm(1)
        X[j+1] <- X[j] + Dt*mu[A[j+1]]*X[j] + sigma[A[j+1]]*X[j]*Winc
    }
    return(list('X' = X, 'A' = A))
}
