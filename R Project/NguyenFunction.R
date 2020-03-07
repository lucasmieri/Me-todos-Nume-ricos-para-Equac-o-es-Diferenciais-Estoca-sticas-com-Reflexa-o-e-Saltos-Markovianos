runNguyen <- function(Dt,L)
{
    # Initializing
    X   = rep(0,floor(L)+1)
    X[1]= X0
    A   = rep(0,floor(L)+1)
    A[1] = sample(c(1,2),1,prob=p0)

    time = seq(0,T,Dt)
    tau <- rexp(1,lambda[A[1]])
        
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
        
        Winc   = sqrt(Dt)*rnorm(1)
        X[j+1] = X[j] + Dt*X[j]*mu[A[j+1]] + X[j]*Winc*sigma[A[j+1]]
                 + 0.5*(sigma[A[j+1]])^2*X[j]*(Winc^2- Dt)
    }
    return(list('X' = X, 'A' = A))
}
