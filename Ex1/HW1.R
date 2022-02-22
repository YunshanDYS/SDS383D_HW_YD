library(matlib) # inv() function
library(matrixcalc) # LU decomposition
library(microbenchmark) 

# Set dimensions
p_vec <- c(5,50,100,900)
n_vec <- c(10,60,200,1000)

for (k in 1:3){
  p <- p_vec[k]
  n <- n_vec[k]
  
  # Generate data
  beta <- rnorm(p,0,1)
  X <- matrix(0,nrow=n,ncol=p)
  W <- diag(1,nrow=n,ncol=n)
  y <- matrix(0,nrow=n,ncol=1)
  for (i in 1:n){
    X[i,] <- rnorm(p,0,1)
    y[i] <- sum(X[i,]*beta) + rnorm(1,0,1) #1/sigma_i^2=1
  }
  
  # function with solve()
  solvelm <- function(X,W,y){
    A <- t(X)%*%W%*%X
    b <- t(X)%*%W%*%y
    beta_hat <- solve(A)%*%b
    return(beta_hat)
  }
  
  # # function with inv()
  # invlm <- function(X,W,y){
  #   A <- t(X)%*%W%*%X
  #   b <- t(X)%*%W%*%y
  #   beta_hat <- inv(A)%*%b
  #   return(beta_hat)
  # }
  
  # function with inv()
  invlm <- function(X,W,y){
    beta_hat <- inv(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
    return(beta_hat)
  }
  
  # function with LU decomposition
  LUlm <- function(X,W,y){
    A <- t(X)%*%W%*%X
    b <- t(X)%*%W%*%y
    # LU decomposition
    decomp <- lu.decomposition(A)
    L <- decomp$L
    U <- decomp$U
    y <- forwardsolve(L, b)
    beta_hat <- backsolve(U, y)
    return(beta_hat)
  }
  
  assign(paste0("out",k),microbenchmark(solvelm(X,W,y), invlm(X,W,y), LUlm(X,W,y), times=5)) 
  cat(k,"\t")
}

out1
out2
out3
out4






