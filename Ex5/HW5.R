library(mvtnorm)
library(invgamma)

# Read data
data <- read.csv("cheese.csv", header = T)

store_names <- unique(data$store)
n <- length(store_names)
# the index of the observations in each store
ind <- list()
Ni <- rep(0,n)
for (i in 1:n){
  ind[[i]] <- which(data$store==store_names[i])
  Ni[i] <- length(ind[[i]])
}
ntotal <- nrow(data)

y <- log(data$vol)
X <- cbind(rep(1,ntotal), log(data$price), data$disp, log(data$price)*data$disp)
p <- 4

# MCMC
# Initialization
M <- 5000 # num of iterations
betai <- list()
for (i in 1:n){
  betai[[i]] <- matrix(1, nrow=M, ncol=p) # each beta_i in each row
}
beta <- matrix(1, nrow=M, ncol=p)
sig2 <- matrix(1, nrow=M, ncol=n)
s2 <- matrix(1, nrow=M, ncol=p)
a <- b <- rep(1,M)
Xi_list <- Yi_list <- list()
for (i in 1:n){
  Xi_list[[i]] <- X[ind[[i]],]
  Yi_list[[i]] <- matrix(y[ind[[i]]],ncol=1)
}
for (m in 2:M){
  # update beta_i
  for (i in 1:n){
    Xi <- Xi_list[[i]]
    Yi <- Yi_list[[i]]
    tilde_Sig <- solve(t(Xi)%*%Xi/sig2[(m-1),i] + diag(1/s2[(m-1),]))
    tilde_mu <- tilde_Sig%*%( t(Xi)%*%Yi/sig2[(m-1),i]  + matrix(1/s2[(m-1),]*beta[(m-1),],ncol=1) )
    betai[[i]][m,] <- rmvnorm(1, tilde_mu, tilde_Sig)
  }
  
  # update beta
  bar_beta <- rep(0,p)
  for (i in 1:n){
    bar_beta <- bar_beta + betai[[i]][m,]
  }
  bar_beta <- bar_beta/n
  for (j in 1:p){
    beta[m,j] <- rnorm(1, bar_beta[j], sqrt(sig2/n))
  }
  
  # update sig2
  for (i in 1:n){
    alpha1 <- (a[(m-1)] + Ni[i])/2
    Xi <- Xi_list[[i]]
    Yi <- Yi_list[[i]]
    beta1 <- b[(m-1)]/2 + sum((Yi - Xi%*%betai[[i]][m,])^2)/2
    sig2[m,i] <- rinvgamma(1, alpha1, beta1)
  }
  
  # update s2
  alpha2 <- (n+1)/2
  beta2 <- sapply(betai,FUN = function(x) (x[m,] - beta[m,])^2)
  beta2 <- apply(beta2,1,sum) + rep(1/2,p)
  for(j in 1:p){
    s2[m,j] <- rinvgamma(1, alpha2, beta2[j])
  }
  
  # update a
  a_grid <- seq(0.5,10,length.out=500)
  pa <- sapply(a_grid,FUN= function(x) sum((-x/2-1)*log(sig2[m,])) ) 
  pa <- pa - mean(pa) # to avoid things in exp() going to Inf
  pa <- exp(pa)
  a[m] <- sample(a_grid,1, prob = pa)
  
  # update b
  b_grid <- seq(0.1,10,length.out=500)
  pb <- sapply(b_grid,FUN= function(x) exp(-sum(x/sig2[m,]/2)) ) 
  b[m] <- sample(b_grid,1, prob = pb)
  
  if (m %% 1000 == 0){
    cat(m,"\t")
  }
}

burnin <- 1000
iters <- burnin:M
# trace plots
png("figures/tra_betai.png",width=600,height=600)
plot(iters,betai[[1]][iters,1],xlab="iterations",ylab="first element of beta_i for store 1",bty="l",type="l")
dev.off()
png("figures/tra_beta.png",width=600,height=600)
plot(iters,beta[iters,1],xlab="iterations",ylab="first element of beta",bty="l",type="l")
dev.off()
png("figures/tra_sig2.png",width=600,height=600)
plot(iters,sig2[iters,1],xlab="iterations",ylab="sigma_i^2 for store 1",bty="l",type="l")
dev.off()
png("figures/tra_s2.png",width=600,height=600)
plot(iters,s2[iters,1],xlab="iterations",ylab="s_j^2 for the first covariate",bty="l",type="l")
dev.off()
# plot(iters,a[iters],xlab="iterations",ylab="a",bty="l",type="l")
# plot(iters,b[iters],xlab="iterations",ylab="b",bty="l",type="l")

post_betai <- matrix(0,nrow=n,ncol=p)
for (i in 1:n){
  post_betai[i,] <- apply(betai[[i]][iters,],2,mean)
}

# plot histogram of beta_i
png("figures/hist_betai.png",width=800,height=800)
par(mfrow = c(2,2))
hist(post_betai[,1],xlab = "", main= "Histogram of beta_0i",freq=F)
hist(post_betai[,2],xlab = "", main= "Histogram of beta_1i",freq=F)
hist(post_betai[,3],xlab = "", main= "Histogram of beta_2i",freq=F)
hist(post_betai[,4],xlab = "", main= "Histogram of beta_3i",freq=F)
dev.off()

# scatter plot of sigma^2_i
sig2_post <- apply(sig2[iters,],2,mean)
png("figures/sig.png",width=800,height=800)
plot(x=1:n,y=sig2_post,xlab="stores",ylab="sigma_i^2")
dev.off()

# plot y against estimated y
y_est <- rep(0,ntotal)
for (i in 1:ntotal){
  index <- which(store_names== data$store[i])
  y_est[i] <- sum(post_betai[index,]*X[i,])
}
png("figures/esty.png",width=800,height=800)
plot(x=y,y=y_est,xlab="y",ylab="estimated y")
dev.off()



