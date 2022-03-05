data <- read.csv("greenbuildings.csv",header=TRUE)
# data <- data[which(data$leasing_rate!=0),]
y <- data$Rent*data$leasing_rate/100
X <- cbind(data$green_rating,data$City_Market_Rent,data$age,data$class_a,data$class_b)
# X <- apply(X,2,scale)
X <- cbind(rep(1,nrow(data)),X)
colnames(X) <- c("intercept","green_rating","City_Market_Rent","age","class_a","class_b")
X_raw <- X[,2:ncol(X)]
# y <- scale(y)
y <- matrix(y,ncol=1)

p <- ncol(X)
n <- nrow(X)
# Prior
K_vec <- rep(0.001,p)
K <- diag(K_vec)
m <- rep(1,p)
d <- 1
eta <- 1
#
Km <- matrix(K_vec*m,ncol=1)
Vbeta_inv <- t(X)%*%X + K
Vbeta <- solve(Vbeta_inv)
mubeta <- Vbeta%*%(t(X)%*%y + Km)
dw <- d + n
etaw <- eta + sum(y^2) + sum(m^2*K_vec) - t(mubeta)%*%Vbeta_inv%*%mubeta
## Bayesian credible interval for green_rating
tq <- qt(1-0.025,df=dw)
mu2 <- mubeta[2]
var2 <- Vbeta[2,2]*etaw/dw
cat("credible interval: (",mu2-tq*sqrt(var2),",",mu2+tq*sqrt(var2),")")

fit <- lm(y~X_raw)
confint(fit, 'X_rawgreen_rating', level=0.95)
summary(fit)
mubeta

## residual
res <- y - X%*%mubeta
jpeg("res.jpeg")
hist(res,breaks=50,freq=FALSE)
dev.off()
plot(data$City_Market_Rent,res)

## Gibbs
library(MASS)
M <- 5000 # num of iterations

# inital
h <- 1
beta_samp <- matrix(NA,ncol=p,nrow=M)
w_samp <- rep(NA,M)
lambda_samp <- matrix(NA,ncol=n,nrow=M)
beta <- beta_samp[1,] <- rep(1,p)
w <- w_samp[1] <- 1
lambda <- lambda_samp[1,] <- rep(1,n)
for (m in 2:M){
  lam_X <- sqrt(lambda)*X
  Vbeta_inv <- t(lam_X)%*%lam_X + K
  Vbeta <- solve(Vbeta_inv)
  mubeta <- Vbeta%*%(t(lambda*X)%*%y + Km)
  beta <- beta_samp[m,] <- mvrnorm(1,mubeta,Vbeta/w)
  
  d_w <- n+d
  eta_w <- eta + sum(y^2*lambda) + sum(m^2*K_vec) - t(mubeta)%*%Vbeta_inv%*%mubeta
  w <- w_samp[m] <- rgamma(1,d_w/2,eta_w/2)
  
  b_lam <- (h + w*(y-X%*%beta)^2)/2
  lambda <- lambda_samp[m,] <- rgamma(n,(h+1)/2,b_lam)
  lambda_samp[m,] <- lambda
}

burnin <- 1000
beta_gibbs <- apply(beta_samp[burnin:M,],2,mean)
res_gibbs <- y - X%*%beta_gibbs
jpeg("res_gibbs.jpeg")
hist(res,breaks=50,col=rgb(0,0,1,1/4),freq=FALSE)
hist(res_gibbs,breaks=50,freq=FALSE,col=rgb(1,0,0,1/4),add=T)
dev.off()

beta_samp[burnin:M,]

stud <- res_gibbs/sqrt(apply((w_samp*lambda_samp)[burnin:M,],2,mean))
jpeg("stud.jpeg")
hist(stud,breaks=30,freq=FALSE)
dev.off()

plot(data$City_Market_Rent,lamb_inv)
lamb_inv <- 1/apply(lambda_samp[burnin:M,],2,mean)
jpeg("lamb_inv.jpeg")
hist(lamb_inv,breaks=30,freq=FALSE)
dev.off()



