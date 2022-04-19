## library
library(mvtnorm)
library(reshape2)
library(ggplot2)

## Function of Squared exponential covariance function
Exp2Sigma <- function(x, b, tau1sq, tau2sq){
  eucDist = as.matrix(dist(x,diag=T,upper=T))
  kron.delta = diag(nrow=length(x))
  tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Function of Matern covariance function with par 5/2
M52Sigma <- function(x, b, tau1sq, tau2sq){
  eucDist = as.matrix(dist(x,diag=T,upper=T))
  kron.delta = diag(nrow=length(x))
  tau1sq*(1 + sqrt(5)*eucDist/b + 5*eucDist^2/3/b^2)*exp(-sqrt(5)*eucDist/b) + tau2sq*kron.delta
}

## Function of Gaussian Process
## with mean 0 and squared exponential covariance function
GP <- function(x, b, tau1sq, tau2sq, covf = "Exp2Sigma"){
  if (covf == "Exp2Sigma"){
    Sigma <- Exp2Sigma(x, b, tau1sq, tau2sq)
  }else if (covf == "M52Sigma"){
    Sigma <- M52Sigma(x, b, tau1sq, tau2sq)
  }
  mean <- rep(0, length(x))
  rmvnorm(1, mean, Sigma)
}

## Simulate GP with different b, tau1sq, tau2sq, Exp2Sigma
set.seed(1234)
x <- runif(100)
tau2sq <- 1e-6
tau1sq <- 1
b_seq <- c(0.1, 0.25, 0.5, 0.7, 1)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b_seq[i], tau1sq, tau2sq)
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("b = ", 5),b_seq)) +
  ggtitle("tau1sq = 1, tau2sq = 1e-6") +
  ylim(-3,3)
ggsave("figures/GP_bseq.png", width = 8, height = 6)

tau2sq <- 1e-6
b <- 0.5
tau1sq_seq <- c(0.1, 1, 2, 3, 5)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b, tau1sq_seq[i], tau2sq)
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("tau1sq = ", 5),tau1sq_seq)) +
  ggtitle("b = 0.5, tau2sq = 1e-6") +
  ylim(-3,3)
ggsave("figures/GP_tau1sq_seq.png", width = 8, height = 6)


tau1sq <- 1
b <- 1
tau2sq_seq <- c(0, 0.0001, 0.001, 0.01, 0.1)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b, tau1sq, tau2sq_seq[i])
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("tau2sq = ", 5),tau2sq_seq)) +
  ggtitle("b = 1, tau1sq = 1") +
  ylim(-3,3)
ggsave("figures/GP_tau2sq_seq.png", width = 8, height = 6)


## Simulate GP with different b, tau1sq, tau2sq, M52Sigma
set.seed(1234)
tau2sq <- 1e-6
tau1sq <- 1
b_seq <- c(0.1, 0.25, 0.5, 0.7, 1)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b_seq[i], tau1sq, tau2sq, covf = "M52Sigma")
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("b = ", 5),b_seq)) +
  ggtitle("tau1sq = 1, tau2sq = 1e-6") +
  ylim(-3,3)
ggsave("figures/GP_bseq_M52.png", width = 8, height = 6)

tau2sq <- 1e-6
b <- 0.5
tau1sq_seq <- c(0.1, 1, 2, 3, 5)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b, tau1sq_seq[i], tau2sq, covf = "M52Sigma")
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("tau1sq = ", 5),tau1sq_seq)) +
  ggtitle("b = 0.5, tau2sq = 1e-6") +
  ylim(-3,3)
ggsave("figures/GP_tau1sq_seq_M52.png", width = 8, height = 6)


tau1sq <- 1
b <- 1
tau2sq_seq <- c(0, 0.0001, 0.001, 0.01, 0.1)
y <- matrix(NA, ncol = 5, nrow = 100)
for (i in 1:5){
  y[ ,i] <- GP(x, b, tau1sq, tau2sq_seq[i], covf = "M52Sigma")
}
df <- cbind(x,y)
df <- as.data.frame(df)
df <- melt(df, id = "x")
ggplot(df,            
       aes(x = x,
           y = value,
           color = variable)) +  
  geom_line() +
  scale_color_hue(labels = paste0(rep("tau2sq = ", 5),tau2sq_seq)) +
  ggtitle("b = 1, tau1sq = 1") +
  ylim(-3,3)
ggsave("figures/GP_tau2sq_seq_M52.png", width = 8, height = 6)

### utilities data
data <- read.csv("utilities.csv",header=T)
y <- data$gasbill/data$billingdays
x <- data$temp
n <- length(y)

tau2sq <- 1e-6
tau1sq <- 1
b <- 15
sigma2 <- 0.61
## Function to calculate C(x,x^*)
C_tilde <- function(x, xstar, b, tau1sq, tau2sq){
  eucDistsq = (x - xstar)^2
  tau1sq*exp(-.5*eucDistsq/b^2) + tau2sq*(x == xstar)
}

C <- Exp2Sigma(x, b, tau1sq, tau2sq)
Cinv <- solve(C + sigma2*diag(n))
xstar_seq <- seq(10,80,length.out = 500)
fhat <- var <- rep(0,500)
for (i in 1:500){
  xstar <- xstar_seq[i]
  Ct <- C_tilde(x, xstar, b, tau1sq, tau2sq)
  Cstar <- C_tilde(xstar, xstar, b, tau1sq, tau2sq)
  fhat[i] <- Ct%*%Cinv%*%as.matrix(y,ncol=1)
  var[i] <- Cstar - Ct%*%Cinv%*%matrix(Ct,ncol=1)
}

## make plots
png("figures/CI_GP.png",width=600,height=600)
plot(x,y,ylim=c(-1,10))
lines(xstar_seq,fhat,col="red")
for (i in 1:500){
  segments(xstar_seq[i], fhat[i]-1.96*sqrt(var[i]), xstar_seq[i], fhat[i]+1.96*sqrt(var[i]),col=alpha("red",0.4))
}
dev.off()

## Function to calculate the log marginal distribution
logmar <- function(y, x, b, tau1sq, tau2sq){
  C <- Exp2Sigma(x, b, tau1sq, tau2sq)
  y_mat <- matrix(y,ncol=1)
  lp <- -n/2*log(2*pi) - 1/2*log(det(C + sigma2*diag(n))) - 1/2*t(y_mat)%*%solve(C + sigma2*diag(n))%*%y_mat
  return(lp)
}

## find optimal parameters
tau1sq_seq <- seq(10,80,length.out = 100)
b_seq <- seq(10,80,length.out = 100)
grid <- expand.grid(tau1sq_seq,b_seq)
lp <- rep(NA,nrow(grid))
for (i in 1:nrow(grid)){
  lp[i] <- logmar(y, x, b = grid[i,2], tau1sq = grid[i,1], tau2sq = 1e-6)
}

# plot heatmap
df <- data.frame(tau1sq = grid[,1], b = grid[,2], lp = lp)
ggplot(df, aes(tau1sq, b, fill= lp)) + 
  geom_tile()
ggsave("figures/heatmap.png",width=8,height=6)
df[which.max(lp),]
tau1sq_hat <- df[which.max(lp),1]
b_hat <- df[which.max(lp),2]


## Use the optimal parameters
tau2sq <- 1e-6
tau1sq <- tau1sq_hat
b <- b_hat
sigma2 <- 0.61
## Function to calculate C(x,x^*)
C_tilde <- function(x, xstar, b, tau1sq, tau2sq){
  eucDistsq = (x - xstar)^2
  tau1sq*exp(-.5*eucDistsq/b^2) + tau2sq*(x == xstar)
}

C <- Exp2Sigma(x, b, tau1sq, tau2sq)
Cinv <- solve(C + sigma2*diag(n))
xstar_seq <- seq(10,80,length.out = 500)
fhat <- var <- rep(0,500)
for (i in 1:500){
  xstar <- xstar_seq[i]
  Ct <- C_tilde(x, xstar, b, tau1sq, tau2sq)
  Cstar <- C_tilde(xstar, xstar, b, tau1sq, tau2sq)
  fhat[i] <- Ct%*%Cinv%*%as.matrix(y,ncol=1)
  var[i] <- Cstar - Ct%*%Cinv%*%matrix(Ct,ncol=1)
}

## make plots
png("figures/CI_GP_2.png",width=600,height=600)
plot(x,y,ylim=c(-1,10))
lines(xstar_seq,fhat,col="red")
for (i in 1:500){
  segments(xstar_seq[i], fhat[i]-1.96*sqrt(var[i]), xstar_seq[i], fhat[i]+1.96*sqrt(var[i]),col=alpha("red",0.4))
}
dev.off()

