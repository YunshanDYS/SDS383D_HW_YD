library(reshape2) 
library(ggpubr)
library(ggplot2)

####
## Kernel smoother

## Kernel Function
# s: distance
# h: band width
Kern <- function(s, h){
  K <- dnorm(s/h,0,1)
  return(K)
}

## Weight Function
# x: previous x
# xstar: new x
# h: band width
# Kern: kernel function
weight <- function(x, xstar, h, Kern){
  w <- Kern(x - xstar, h)/h
  w <- w/sum(w)
  return(w)
}

## Kernel Regression Function
Kern_reg <- function(y, x, xstar, h){
  w <- weight(x, xstar, h, Kern)
  fhat <- sum(w*y)
  return(fhat)
}

## Simulation
## y = x^2 + \eps

## Generate data
set.seed(1999)
x <- rnorm(50,0,2)
y <- x^2 + rnorm(50,0,1)

ngrid <- 500
x_grid <- seq(-3,3,length.out = ngrid)
h_array <- c(0.1,0.5,1,2)
f_grid <- matrix(NA,ncol=length(h_array), nrow=ngrid)
for (k in 1:length(h_array)){
  for (i in 1:ngrid){
    f_grid[i,k] <- Kern_reg(y,x,x_grid[i],h_array[k])
  }
}

df <- data.frame(cbind(x_grid, x_grid^2, f_grid))
colnames(df) <- c("x", "f(x)", paste0(rep("h = ",length(h_array)), h_array))
df_plot <- melt(df, id = "x")

ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line()
ggsave("figures/kernel.png", width=8, height=6)

### Cross Validation
nh <- 100
h_seq <- seq(0.1,2,length.out=nh)
ngrid <- 500
x_grid <- seq(-3,3,length.out = ngrid)

## Cross Validation Function
cross_va <- function(x_train, y_train, x_test, y_test){
  f_grid <- matrix(NA,ncol=nh, nrow=ngrid)
  y_est <- matrix(NA,ncol=nh, nrow=length(y_test))
  err <- rep(NA,nh)
  for (k in 1:nh){
    for (i in 1:ngrid){
      f_grid[i,k] <- Kern_reg(y_train,x_train,x_grid[i],h_seq[k])
    }
    for (i in 1:length(y_test)){
      y_est[i,k] <- Kern_reg(y_train,x_train,x_test[i],h_seq[k])
    }
    err[k] <- mean((y_est[,k] - y_test)^2)
  }
  return(list(f_grid = f_grid, err = err))
}

## smooth, low noise
## f(x) = x^2, eps \sim N(0,0.5^2)
x <- rnorm(1000,0,2)
y <- x^2 + rnorm(1000,0,0.5)
out1 <- cross_va(x_train=x[1:500], y_train=y[1:500], x_test=x[501:1000], y_test=y[501:1000])
ind <- which.min(out1$err)
h_sel <- h_seq[ind]
f_grid <- out1$f_grid[,ind]
df <- data.frame(cbind(x_grid, x_grid^2, f_grid, out1$f_grid[,c(25,50,75)]))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)), paste0(rep("h = ",3), round(h_seq[c(25,50,75)],4)))
df_plot <- melt(df, id = "x")

p1 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Smooth, low noise")
p1

## smooth, high noise
## f(x) = x^2, eps \sim N(0,5^2)
x <- rnorm(1000,0,2)
y <- x^2 + rnorm(1000,0,2)
out1 <- cross_va(x_train=x[1:500], y_train=y[1:500], x_test=x[501:1000], y_test=y[501:1000])
ind <- which.min(out1$err)
h_sel <- h_seq[ind]
f_grid <- out1$f_grid[,ind]
df <- data.frame(cbind(x_grid, x_grid^2, f_grid, out1$f_grid[,c(25,50,75)]))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)), paste0(rep("h = ",3), round(h_seq[c(25,50,75)],4)))
df_plot <- melt(df, id = "x")

p2 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Smooth, high noise")
p2

## wiggly, low noise
x <- rnorm(1000,0,2)
y <- sin(5*x) + rnorm(1000,0,0.1)
out1 <- cross_va(x_train=x[1:500], y_train=y[1:500], x_test=x[501:1000], y_test=y[501:1000])
ind <- which.min(out1$err)
h_sel <- h_seq[ind]
f_grid <- out1$f_grid[,ind]
df <- data.frame(cbind(x_grid, sin(5*x_grid), f_grid, out1$f_grid[,c(25,50,75)]))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)), paste0(rep("h = ",3), round(h_seq[c(25,50,75)],4)))
df_plot <- melt(df, id = "x")

p3 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Wiggly, low noise")
p3

## wiggly, high noise
x <- rnorm(1000,0,2)
y <- sin(5*x) + rnorm(1000,0,1)
out1 <- cross_va(x_train=x[1:500], y_train=y[1:500], x_test=x[501:1000], y_test=y[501:1000])
ind <- which.min(out1$err)
h_sel <- h_seq[ind]
f_grid <- out1$f_grid[,ind]
df <- data.frame(cbind(x_grid, sin(5*x_grid), f_grid, out1$f_grid[,c(25,50,75)]))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)), paste0(rep("h = ",3), round(h_seq[c(25,50,75)],4)))
df_plot <- melt(df, id = "x")

p4 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Wiggly, high noise")
p4

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggsave("figures/crossva.png",width=15,height=12)

## LOOCV
## Cross Validation Function
nh <- 100
h_seq <- seq(0.1,2,length.out=nh)
ngrid <- 500
x_grid <- seq(-3,3,length.out = ngrid)
cross_va_LOO <- function(x, y){
  f_grid <- rep(NA,ngrid)
  err <- rep(NA,nh)
  n <- length(y)
  H <- matrix(NA,ncol=n, nrow=n)
  for (k in 1:nh){
    for (i in 1:n){
      H[i,] <- weight(x, x[i], h_seq[k], Kern)
    }
    yhat <- H%*%y
    err[k] <- sum((y-yhat/(1-diag(H)))^2)
  }
  h <- h_seq[which.min(err)]
  for (i in 1:ngrid){
    f_grid[i] <- Kern_reg(y,x,x_grid[i],h)
  }
  return(list(f_grid = f_grid, h = h))
}

## smooth, low noise
## f(x) = x^2, eps \sim N(0,0.5^2)
x <- rnorm(1000,0,2)
y <- x^2 + rnorm(1000,0,0.5)
out1 <- cross_va_LOO(x,y)
h_sel <- out1$h
f_grid <- out1$f_grid
df <- data.frame(cbind(x_grid, x_grid^2, f_grid))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)))
df_plot <- melt(df, id = "x")

p1 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Smooth, low noise")
p1

## smooth, high noise
## f(x) = x^2, eps \sim N(0,5^2)
x <- rnorm(1000,0,2)
y <- x^2 + rnorm(1000,0,2)
out1 <- cross_va_LOO(x,y)
h_sel <- out1$h
f_grid <- out1$f_grid
df <- data.frame(cbind(x_grid, x_grid^2, f_grid))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)))
df_plot <- melt(df, id = "x")

p2 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Smooth, high noise")
p2

## wiggly, low noise
x <- rnorm(1000,0,2)
y <- sin(5*x) + rnorm(1000,0,0.1)
out1 <- cross_va_LOO(x,y)
h_sel <- out1$h
f_grid <- out1$f_grid
df <- data.frame(cbind(x_grid, sin(5*x_grid), f_grid))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)))
df_plot <- melt(df, id = "x")

p3 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Wiggly, low noise")
p3

## wiggly, high noise
x <- rnorm(1000,0,2)
y <- sin(5*x) + rnorm(1000,0,1)
out1 <- cross_va_LOO(x,y)
h_sel <- out1$h
f_grid <- out1$f_grid
df <- data.frame(cbind(x_grid, sin(5*x_grid), f_grid))
colnames(df) <- c("x", "f(x)", paste0("Selected h = ", round(h_sel,4)))
df_plot <- melt(df, id = "x")

p4 <- ggplot(df_plot, aes(x = x, y = value, color = variable)) +
  geom_line() +
  ggtitle("Wiggly, high noise")
p4

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggsave("figures/crossva_LOO.png",width=15,height=12)


### Local linear regression
weight <- function(x, xstar, h, Kern){
  w <- Kern(x - xstar, h)
  return(w)
}

## Function for s_1(x)
s1 <- function(x, x_new, h){
  n <- length(x_new)
  s1 <- rep(0,n)
  for (i in 1:n){
    s1[i] <- sum(weight(x, x_new[i], h, Kern) * (x - x_new[i]))
  }
  return(s1)
}

## Function for s_2(x)
s2 <- function(x, x_new, h){
  n <- length(x_new)
  s2 <- rep(0,n)
  for (i in 1:n){
    s2[i] <- sum(weight(x, x_new[i], h, Kern) * (x - x_new[i])^2)
  }
  return(s2)
}

## Function for w_i(x) (H matrix)
wi <- function(x, x_new, h){
  n <- length(x_new)
  wi <- matrix(0,ncol=n,nrow=n)
  for (i in 1:n){
    wi[i,] <- weight(x, x_new[i], h, Kern) * (s2(x,x_new,h) - (x - x_new[i])*s1(x,x_new,h))
    wi[i,] <- wi[i,]/sum(wi[i,])
  }
  return(wi)
}

## Cross Validation Function
nh <- 100
h_seq <- seq(0.1,10,length.out=nh)
cross_va_LOO_2 <- function(x, y){
  err <- rep(NA,nh)
  n <- length(y)
  for (k in 1:nh){
    H <- wi(x, x, h_seq[k])
    yhat <- H%*%y
    err[k] <- sum((y-yhat/(1-diag(H)))^2)
  }
  h <- h_seq[which.min(err)]
  H <- wi(x, x, h)
  yhat <- H%*%y
  return(list(yhat = yhat, h = h, H=H))
}

data <- read.csv("utilities.csv",header=T)
y <- data$gasbill/data$billingdays
x <- data$temp

# x <- scale(x)
# y <- scale(y)

out <- cross_va_LOO_2(x,y)

yhat <- out$yhat
out$h
H <- out$H
png("figures/CI.png",width=600,height=600)
plot(x,y)
points(x,yhat,col="red")
sig2 <- as.numeric(var(y))
n <- length(y)
sd <- rep(0,n)
for (i in 1:n){
  sd[i] <- sqrt(sig2*sum(H[i,]^2))
  segments(x[i], yhat[i]-2*sd[i], x[i], yhat[i]+2*sd[i],col="red")
}
dev.off()
## residuals
png("figures/residuals.png",width=600,height=600)
residuals <- yhat-y
plot(x,residuals)
dev.off()





