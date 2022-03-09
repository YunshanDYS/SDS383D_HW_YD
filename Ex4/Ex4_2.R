library(dplyr)
library(invgamma)
data <- read.csv("bloodpressure.csv",header=TRUE)

y1 <- data$systolic[which(data$treatment==1)]
y2 <- data$systolic[which(data$treatment==2)]
test1 <- t.test(y1,y2)
test1
test1$stderr

bar_y <- data %>% group_by(subject,treatment) %>% summarize(ybar=mean(systolic))
y1 <- bar_y$ybar[which(bar_y$treatment==1)]
y2 <- bar_y$ybar[which(bar_y$treatment==2)]
test2 <- t.test(y1,y2)
test2
test2$stderr

######
y_subj <- data %>%
  group_by(subject) %>%
  summarize(bar_y = mean(systolic), n_data = n(), xi = as.numeric(treatment[1]==2))
y_subj

P <- nrow(y_subj)
Ni <- y_subj$n_data
ybar <- y_subj$bar_y
xi <- y_subj$xi

## Gibbs
set.seed(1234)
M <- 10000
# initialize
theta <- matrix(NA,ncol=P,nrow=M)
mu <- sigma2 <- tau2 <- beta <- rep(NA,M)
beta[1] <- 1
theta[1,] <- rep(1,P)
mu[1] <- sigma2[1] <- tau2[1] <- 1

for (m in 2:M){
  # update beta
  sig2_beta <- 1/(sum(xi^2)/tau2[m-1]/sigma2[m-1] + 1)
  mu_beta <- sig2_beta*(sum((theta[m-1,]-mu[m-1])*xi)/tau2[m-1]/sigma2[m-1])
  beta[m] <- rnorm(1,mean=mu_beta,sd=sqrt(sig2_beta))
  
  # update theta
  tilde_sig <- 1/(Ni/sigma2[m-1] + 1/(tau2[m-1]*sigma2[m-1]))
  tilde_mu <- tilde_sig*(Ni*ybar/sigma2[m-1] + (mu[m-1] + beta[m]*xi)/(tau2[m-1]*sigma2[m-1]))
  for (i in 1:P){
    theta[m,i] <- rnorm(1,mean=tilde_mu[i], sd=sqrt(tilde_sig[i]))
  }
  
  # update mu
  mu[m] <- rnorm(1,mean=mean(theta[m,] - beta[m]*xi),sd=sqrt(tau2[m-1]*sigma2[m-1]/P))
  
  # update sigma^2
  b_sig <- 0
  for (i in 1:P){
    b_sig <- b_sig + sum((data[which(data$subject==i),2] - theta[m,i])^2)/2 + (theta[m,i] - mu[m] - beta[m]*xi[i])^2/tau2[m-1]/2
  }
  a_sig <- (sum(Ni) +P)/2
  sigma2[m] <- rinvgamma(1,a_sig,b_sig)
  
  # update tau^2
  b_tau <- 1/2 + sum((theta[m,] - mu[m] - beta[m]*xi)^2)/sigma2[m]/2
  a_tau <- (P+1)/2
  tau2[m] <- rinvgamma(1,a_tau,b_tau)
}

burnin <- 1000
iters <- burnin:M
jpeg("figures/2mu_tra.jpeg",width=600,height=600)
plot(iters,mu[iters],xlab="iterations",ylab=expression(mu),bty="l",type="l")
dev.off()
jpeg("figures/2beta_tra.jpeg",width=600,height=600)
plot(iters,beta[iters],xlab="iterations",ylab=expression(beta),bty="l",type="l")
dev.off()
jpeg("figures/2sig_tra.jpeg",width=600,height=600)
plot(iters,sigma2[iters],xlab="iterations",ylab=expression(sigma^2),bty="l",type="l")
dev.off()
jpeg("figures/2tau_tra.jpeg",width=600,height=600)
plot(iters,tau2[iters],xlab="iterations",ylab=expression(tau^2),bty="l",type="l")
dev.off()

beta_post_mean <- mean(beta[iters])
beta_post_mean
beta_post_sd <- sd(beta[iters])
beta_post_sd



