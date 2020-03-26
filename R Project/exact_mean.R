library(expm)

exact_mean <- function()
{
    A <- matrix(c(mu[1] - lambda[1], lambda[2], lambda[1],mu[2]-lambda[2]),nrow=2)
    C <- c(X0*p0[1],X0*p0[2])

    u1  <- C %*% expm(A*1)
    EX1 <- u1[1] + u1[2]
    return(EX1)
}