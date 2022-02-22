## Gradient descent

# Data
data <- read.csv("wdbc.csv",header=FALSE)

y <- as.numeric(data[,2]=="M")
X <- as.matrix(data[,3:12])
X <- apply(X,2,scale)
X <- cbind(rep(1,nrow(data)),X)

p <- 11

# log likelihood
loglike <- function(y,theta){
  ll <-  sum(y*theta - log(1+exp(theta)))
  return(ll)
}

# Evaluates the gradient of the negative
# log likelihood of the logit model
# @beta: current guess for coefficients, including intercept
# @y: successes
# @m: sample sizes
# @X: feature matrix (with column of 1's for intercept)
logit_nll_grad = function(beta, y, X, m) {
  psi_hat = drop(X %*% beta)
  y_hat = m*1/{1+exp(-psi_hat)}
  delta = y_hat - y
  -crossprod(X, delta)
}

## Fixed step size
# Initialization
beta_iter <- beta <- matrix(rep(0.1,11),ncol=1)
ll_iter <- c()
tol <- 1e-6
gamma <- 0.01
conv <- FALSE
niter <- 0
ll_old <- -1000
while (!conv){
  g <- logit_nll_grad(beta, y, X, m = rep(1, nrow(X)))
  beta <- beta + gamma*g
  beta_iter <- cbind(beta_iter,beta)
  theta <- X%*%beta
  ll <- loglike(y,theta)
  ll_iter <- c(ll_iter,ll)
  if (abs((ll_old - ll)/ll_old) < tol){
    conv <- TRUE
  }
  niter <- niter + 1
  ll_old <- ll
}

ll_iter0 <- ll_iter

# Line search
f_gamma <- function(gamma,beta,g,y,X){
  theta <- X%*%(beta+gamma*g)
  ll <-  sum(y*theta - log(1+exp(theta)))
  return(-ll)
}

beta_iter <- beta <- matrix(rep(0.1,11),ncol=1)
ll_iter <- c()
tol <- 1e-6
conv <- FALSE
niter <- 0
ll_old <- -1000
while (!conv){
  g <- logit_nll_grad(beta, y, X, m = rep(1, nrow(X)))
  gamma <- optimize(f_gamma,c(1e-6,0.3),beta=beta,g=g,y=y,X=X)$minimum
  beta <- beta + gamma*g
  beta_iter <- cbind(beta_iter,beta)
  theta <- X%*%beta
  ll <- loglike(y,theta)
  ll_iter <- c(ll_iter,ll)
  if (abs((ll_old - ll)/ll_old) < tol){
    conv <- TRUE
  }
  niter <- niter + 1
  ll_old <- ll
}

ll_iter_line <- ll_iter

## Check by comparing to glm()
data_glm <- data.frame(cbind(y,X))

fit <- glm(y~V3+V4+V5+V6+V7+V8+V9+V10+V11+V12,data=data_glm,family=binomial())
est_glm <- fit$coefficients
beta_glm <- matrix(est_glm,ncol=1)
theta_glm <- X%*%beta_glm
ll_glm <- loglike(y,theta_glm)

# plot
jpeg("logl.jpg",width=600,height=600)
plot(ll_iter_line,cex=0.3, ylab="log likelihood",xlab="iterations")
points(ll_iter0,cex=0.3, ylab="log likelihood",xlab="iterations",col="blue")
abline(h = ll_glm, col ="red")
legend("bottomright", legend=c("glm", "line search","fixed step"),
       col=c("red","black", "blue"), lty=1:2, cex=0.8)
dev.off()

### Newton
Heis <- function(beta, y, X){
  psi_hat <- drop(X %*% beta)
  mu <- 1/(exp(-psi_hat)+1)
  # W <- diag(1/mu/(1-mu))
  # H <- t(X)%*%W%*%X
  W <- mu*(1-mu)
  xx <- X*sqrt(W)
  H <- -t(xx)%*%xx
  return(H)
}

beta_iter <- beta <- matrix(rep(0.1,11),ncol=1)
ll_iter <- c()
tol <- 1e-6
conv <- FALSE
niter <- 0
ll_old <- -1000
while (!conv){
  g <- logit_nll_grad(beta, y, X, m = rep(1, nrow(X)))
  H <- Heis(beta,y,X)
  beta <- beta - solve(H)%*%g
  beta_iter <- cbind(beta_iter,beta)
  theta <- X%*%beta
  ll <- loglike(y,theta)
  ll_iter <- c(ll_iter,ll)
  if (abs((ll_old - ll)/ll_old) < tol){
    conv <- TRUE
  }
  niter <- niter + 1
  ll_old <- ll
}

# number of iteration to converge
cat("number of iteration to converge: ",niter)

## Check by comparing to glm()
data_glm <- data.frame(cbind(y,X))

fit <- glm(y~V3+V4+V5+V6+V7+V8+V9+V10+V11+V12,data=data_glm,family=binomial)
est_glm <- fit$coefficients
beta_glm <- matrix(est_glm,ncol=1)
theta_glm <- X%*%beta_glm
ll_glm <- loglike(y,theta_glm)

beta
# difference between glm and newton
beta_glm - beta

# plot
jpeg("logl_2.jpg",width=600,height=600)
plot(ll_iter_line,cex=0.3, ylab="log likelihood",xlab="iterations")
points(ll_iter0,cex=0.3, ylab="log likelihood",xlab="iterations",col="blue")
points(ll_iter,cex=0.3, ylab="log likelihood",xlab="iterations",col="green")
lines(ll_iter,cex=0.3, ylab="log likelihood",xlab="iterations",col="green")
abline(h = ll_glm, col ="red")
legend("bottomright", legend=c("glm", "line search","fixed step","Newton"),
       col=c("red","black", "blue","green"), lty=1:2, cex=0.8)
dev.off()

ll_iter

# sqrt of diag of inverse negative Hessian
sqrt(diag(solve(-H)))
std_glm <- summary(fit)$coefficients[,2]
# std error in glm()
std_glm
std_glm - sqrt(diag(solve(-H)))




