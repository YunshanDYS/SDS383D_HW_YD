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
set.seed(1999)
x <- runif(100)
tau2sq <- 1e-6
tau1sq <- 1
b_seq <- c(0.01, 0.1, 0.5, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("b = ", 5),b_seq)) +
  ggtitle("tau1sq = 1, tau2sq = 1e-6")
ggsave("figures/GP_bseq.png", width = 8, height = 6)

tau2sq <- 1e-6
b <- 1
tau1sq_seq <- c(0.01, 0.1, 0.5, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("tau1sq = ", 5),tau1sq_seq)) +
  ggtitle("b = 1, tau2sq = 1e-6")
ggsave("figures/GP_tau1sq_seq.png", width = 8, height = 6)


tau1sq <- 1
b <- 1
tau2sq_seq <- c(0, 0.01, 0.1, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("tau2sq = ", 5),tau2sq_seq)) +
  ggtitle("b = 1, tau1sq = 1")
ggsave("figures/GP_tau2sq_seq.png", width = 8, height = 6)


## Simulate GP with different b, tau1sq, tau2sq, M52Sigma
set.seed(1999)
tau2sq <- 1e-6
tau1sq <- 1
b_seq <- c(0.01, 0.1, 0.5, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("b = ", 5),b_seq)) +
  ggtitle("tau1sq = 1, tau2sq = 1e-6")
ggsave("figures/GP_bseq_M52.png", width = 8, height = 6)

tau2sq <- 1e-6
b <- 1
tau1sq_seq <- c(0.01, 0.1, 0.5, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("tau1sq = ", 5),tau1sq_seq)) +
  ggtitle("b = 1, tau2sq = 1e-6")
ggsave("figures/GP_tau1sq_seq_M52.png", width = 8, height = 6)


tau1sq <- 1
b <- 1
tau2sq_seq <- c(0, 0.01, 0.1, 1, 2)
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
  geom_point() +
  scale_color_hue(labels = paste0(rep("tau2sq = ", 5),tau2sq_seq)) +
  ggtitle("b = 1, tau1sq = 1")
ggsave("figures/GP_tau2sq_seq_M52.png", width = 8, height = 6)






