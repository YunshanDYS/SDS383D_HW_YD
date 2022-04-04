library(mvtnorm)
# library(invgamma)
library(LaplacesDemon)
library(truncnorm)
library(ggplot2)
library(ggpubr)

# Read data
data <- read.csv("polls.csv", header = T)
data <- data[!is.na(data$bush),4:10]
ntotal <- nrow(data)

# change categorical edu and age to numbers
edu_names <- unique(data$edu)
edu_names
age_names <- unique(data$age)
age_names
data$edu_NoHS <- (data$edu=="NoHS")
data$edu_HS <- (data$edu=="HS")
data$edu_SomeColl <- (data$edu=="SomeColl")
data$edu_Bacc <- (data$edu=="Bacc")
data$age_18to29 <- (data$age=="18to29")
data$age_30to44 <- (data$age=="30to44")
data$age_45to64 <- (data$age=="45to64")
data$age_65plus <- (data$age=="65plus")
# age_names <- sort(age_names)
# for (i in 1:ntotal){
#   data[i,3] <- which(edu_names==data[i,3])
#   data[i,4] <- which(age_names==data[i,4])
# }

state_names <- unique(data$state)
n <- length(state_names)
# the index of the observations in each store
ind <- list()
Ni <- rep(0,n)
for (i in 1:n){
  ind[[i]] <- which(data$state==state_names[i])
  Ni[i] <- length(ind[[i]])
}


y <- as.numeric(data$bush)
X <- cbind(rep(1,ntotal), 
           data$edu_NoHS, data$edu_HS, data$edu_SomeColl, #data$edu_Bacc,
           data$age_18to29, data$age_30to44, data$age_45to64, #data$age_65plus,
           data$female, data$black, scale(data$weight))
# X <- cbind(rep(1,ntotal),
#            data$edu_NoHS, data$edu_HS, data$edu_SomeColl, data$edu_Bacc,
#            data$age_18to29, data$age_30to44, data$age_45to64, data$age_65plus,
#            data$female, data$black, scale(data$weight))
X <- apply(X,2,as.numeric)
p <- ncol(X)

# MCMC
# Initialization
M <- 5000 # num of iterations
betai <- list()
for (i in 1:n){
  betai[[i]] <- matrix(1, nrow=M, ncol=p) # each beta_i in each row
}
aij <- matrix(1, nrow=M, ncol=ntotal)
beta <- matrix(1, nrow=M, ncol=p)
s2 <- matrix(1, nrow=M, ncol=p)
Xi_list <- Yi_list <- list()
for (i in 1:n){
  Xi_list[[i]] <- X[ind[[i]],]
  Yi_list[[i]] <- matrix(y[ind[[i]]],ncol=1)
}
Sigma <- diag(p) ##
nu.0 <- p+1 ##
for (m in 2:M){
  # # update beta_i
  # for (i in 1:n){
  #   Xi <- Xi_list[[i]]
  #   ai <- matrix(aij[(m-1),][ind[[i]]],ncol=1)
  #   tilde_Sig <- solve(t(Xi)%*%Xi + diag(1/s2[(m-1),]))
  #   tilde_mu <- tilde_Sig%*%( t(Xi)%*%ai  + matrix(1/s2[(m-1),]*beta[(m-1),],ncol=1) )
  #   betai[[i]][m,] <- rmvnorm(1, tilde_mu, tilde_Sig)
  # }
  
  for(i in 1:n){
    Xi <- Xi_list[[i]]
    ai <- matrix(aij[(m-1),][ind[[i]]],ncol=1)
    B.tilde <- solve(solve(Sigma) + t(Xi)%*%Xi)
    beta.tilde <- B.tilde%*%(  solve(Sigma)%*%matrix(beta[(m-1),],ncol=1) +  t(Xi)%*%ai )

    betai[[i]][m,] <- rmvnorm(1, beta.tilde, B.tilde, checkSymmetry = FALSE) # added last argument
  }
  
  # update beta
  bar_beta <- rep(0,p)
  for (i in 1:n){
    bar_beta <- bar_beta + betai[[i]][m,]
  }
  bar_beta <- bar_beta/n
  for (j in 1:p){
    beta[m,j] <- rnorm(1, bar_beta[j], sqrt(s2[m-1,j]/n))
  }
  
  # # update s2
  # alpha2 <- (n+1)/2
  # beta2 <- sapply(betai,FUN = function(x) (x[m,] - beta[m,])^2)
  # beta2 <- apply(beta2,1,sum)/2 + rep(1/2,p)
  # for(j in 1:p){
  #   s2[m,j] <- rinvgamma(1, alpha2, beta2[j])
  # }

  ###

  nu <- n + nu.0
  tmp.sum <- 0
  for(i in 1:n){
    tmp.sum <- tmp.sum + (betai[[i]][m,] - beta[m,])%*%t(betai[[i]][m,] - beta[m,])
  }
  S <- diag(p) + tmp.sum

  Sigma <- rinvwishart(nu, S)
  
  # update aij
  for (i in 1:ntotal){
    if (y[i]==1){
      aij[m,i] <- rtruncnorm(1, a = 0, mean = sum(X[i,]*betai[[which(state_names==data$state[i])]][m,]), sd = 1)
      # aij_grid <- seq(1e-4,5,length.out=500)
      # pa <- dnorm(aij_grid, sum(X[i,]*betai[[which(state_names==data$state[i])]][m,]),1)
      # aij[m,i] <- sample(aij_grid,1,prob=pa)
    }else{
      aij[m,i] <- rtruncnorm(1, b = 0, mean = sum(X[i,]*betai[[which(state_names==data$state[i])]][m,]), sd = 1)
      # aij_grid <- seq(-5,-1e-4,length.out=500)
      # pa <- dnorm(aij_grid, sum(X[i,]*betai[[which(state_names==data$state[i])]][m,]),1)
      # aij[m,i] <- sample(aij_grid,1,prob=pa)
    }
  }
  
  
  if (m %% 1000 == 0){
    cat(m,"\t")
  }
}

burnin <- 2500
iters <- burnin:M
# trace plots
png("figures/2tra_betai.png",width=600,height=600)
plot(iters,betai[[1]][iters,1],xlab="iterations",ylab="first element of beta_i for state 1",bty="l",type="l")
dev.off()
png("figures/2tra_beta.png",width=600,height=600)
plot(iters,beta[iters,1],xlab="iterations",ylab="first element of beta",bty="l",type="l")
dev.off()
png("figures/2tra_aij.png",width=600,height=600)
plot(iters,aij[iters,1],xlab="iterations",ylab="a_ij for state 1",bty="l",type="l")
dev.off()
png("figures/2tra_s2.png",width=600,height=600)
plot(iters,s2[iters,1],xlab="iterations",ylab="s_j^2 for the first covariate",bty="l",type="l")
dev.off()

post_betai <- matrix(0,nrow=n,ncol=p)
for (i in 1:n){
  post_betai[i,] <- apply(betai[[i]][iters,],2,mean)
}

# # plot histogram of beta_i
# png("figures/2hist_betai.png",width=800,height=800)
# par(mfrow = c(2,2))
# hist(post_betai[,1],xlab = "", main= "Histogram of beta_0i",freq=F)
# hist(post_betai[,2],xlab = "", main= "Histogram of beta_1i",freq=F)
# hist(post_betai[,3],xlab = "", main= "Histogram of beta_2i",freq=F)
# hist(post_betai[,4],xlab = "", main= "Histogram of beta_3i",freq=F)
# dev.off()

plot_title <- c("intercept", "edu_NoHS", "edu_HS", "edu_SomeColl", 
                "age_18to29", "age_30to44", "age_45to64", "female", "black", "weight")
post_beta <- apply(beta[iters,],2,mean)
df_plot <- data.frame(covariates = plot_title, post = post_beta, col= as.character(as.numeric(post_beta > 0)) )
p0 <- ggplot() +
  geom_point(data = df_plot, aes(x = covariates, y = post, color = col),size=2) +
  scale_color_manual(values = c("1" = "red", "0" = "blue")) +
  theme(legend.position="none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") 
p0
ggsave(paste0("figures/2post covariates0.png"),width=12,height=8)
p0

df_plot <- data.frame(state = state_names, post = post_beta[,i], col= as.character(as.numeric(post_betai[,i] > 0)) )
p0 <- ggplot() +
  geom_point(data = df_plot, aes(x = state, y = post, color = col),size=2) +
  scale_color_manual(values = c("1" = "red", "0" = "blue")) +
  theme(legend.position="none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle(plot_title[i]) 
assign(paste0("p",i),p0)
ggsave(paste0("figures/2post covariates",i, ".png"),width=12,height=8)

# par(mfrow = c(3,3))
# plot(1:n,post_betai[,1],xlab="States",ylab="mu", xlim=state_names, ylim=-2,2)
# plot(1:n,post_betai[,2],xlab="States",ylab="beta1")
# plot(1:n,post_betai[,3],xlab="States",ylab="beta2")
# plot(1:n,post_betai[,4],xlab="States",ylab="beta3")
# plot(1:n,post_betai[,5],xlab="States",ylab="beta4")
# plot(1:n,post_betai[,6],xlab="States",ylab="beta5")
# plot(1:n,post_betai[,7],xlab="States",ylab="beta6")
# plot(1:n,post_betai[,8],xlab="States",ylab="beta7")
# plot(1:n,post_betai[,9],xlab="States",ylab="beta8")
# dev.off()
# 
# state_names[which(post_betai[,1]<0)]

# # scatter plot of sigma^2_i
# sig2_post <- apply(s2[iters,],2,mean)
# png("figures/2sig.png",width=800,height=800)
# plot(x=1:p,y=sig2_post,xlab="stores",ylab="sigma_k^2")
# dev.off()
# 
# # plot y against estimated y
# y_est <- rep(0,ntotal)
# for (i in 1:ntotal){
#   index <- which(state_names== data$state[i])
#   y_est[i] <- (sum(post_betai[index,]*X[i,])>0)
# }
# png("figures/2esty.png",width=800,height=800)
# plot(x=y,y=y_est-y,xlab="y",ylab="estimated y")
# dev.off()



