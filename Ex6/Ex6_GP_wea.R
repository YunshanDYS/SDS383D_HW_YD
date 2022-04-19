####
data <- read.csv("weather.csv",header=T)
y <- data$pressure
# y <- data$temperature
X <- data[,3:4]
n <- length(y)

## Function of Squared exponential covariance function
Exp2Sigma <- function(X, b, tau1sq, tau2sq){
  eucDist = as.matrix(dist(X,diag=T,upper=T))
  kron.delta = diag(nrow=nrow(X))
  tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Function of Squared exponential covariance function
Exp2Sigma_sep <- function(X, b, tau1sq, tau2sq){
  Dist1 = as.matrix(dist(X[,1],diag=T,upper=T))
  Dist2 = as.matrix(dist(X[,2],diag=T,upper=T))
  kron.delta = diag(nrow=nrow(X))
  tau1sq*exp(-.5*((Dist1/as.numeric(b[1]))^2+ (Dist2/as.numeric(b[2]))^2)) + tau2sq*kron.delta
}

## Function to calculate the log marginal distribution
logmar <- function(y, X, b, tau1sq, tau2sq){
  C <- Exp2Sigma(X, b, tau1sq, tau2sq)
  y_mat <- matrix(y,ncol=1)
  lp <- -n/2*log(2*pi) - 1/2*log(det(C + sigma2*diag(n))) - 1/2*t(y_mat)%*%solve(C + sigma2*diag(n))%*%y_mat
  return(lp)
}

tau2sq <- 1e-6
sigma2 <- 1
## find optimal parameters
# tau1sq_seq <- seq(1,20,length.out = 20) ## for temp
tau1sq_seq <- seq(100,300,length.out = 20) ## for pressure
b1_seq <- seq(0.01,5,length.out = 20)
grid <- expand.grid(tau1sq_seq,b1_seq)
lp <- rep(NA,nrow(grid))
for (i in 1:nrow(grid)){
  lp[i] <- logmar(y, X, b = grid[i,2], tau1sq = grid[i,1], tau2sq = 1e-6)
  if (i %%100 == 0){cat(i,"\t")}
}

# plot heatmap
df <- data.frame(tau1sq = grid[,1], b1 = grid[,2], lp = lp)
# ggplot(df, aes(tau1sq, b, fill= lp)) + 
#   geom_tile()
# ggsave("figures/heatmap_wea.png",width=8,height=6)
df[which.max(lp),]
tau1sq_hat <- df[which.max(lp),1]
b_hat <- df[which.max(lp),2]

### Use the optimal parameters
n_grid <- 100
x_long <- seq(min(X[,1]), max(X[,1]), length.out = n_grid)
x_lat <- seq(min(X[,2]), max(X[,2]), length.out = n_grid)
X_new <- expand.grid(x_long, x_lat) 
tau2sq <- 1e-6
tau1sq <- tau1sq_hat
b <- b_hat
sigma2 <- 1
## Function to calculate C(x,x^*)
C_tilde <- function(X, Xstar, b, tau1sq, tau2sq){
  Dist1 = (X[,1] - Xstar[,1])
  Dist2 = (X[,2] - Xstar[,2])
  tau1sq*exp(-.5*((Dist1/as.numeric(b))^2+ (Dist2/as.numeric(b))^2)) #+ tau2sq*colSums(X == Xstar)
}

C <- Exp2Sigma(X, b, tau1sq, tau2sq)
Cinv <- solve(C + sigma2*diag(n))
fhat <- var <- rep(0,nrow(X_new))
for (i in 1:nrow(X_new)){
  Ct <- C_tilde(X, X_new[i,], b, tau1sq, tau2sq)
  Cstar <- C_tilde(X_new[i,], X_new[i,], b, tau1sq, tau2sq)
  fhat[i] <- Ct%*%Cinv%*%as.matrix(y,ncol=1)
  var[i] <- Cstar - Ct%*%Cinv%*%matrix(Ct,ncol=1)
}

df <- data.frame(x=X_new[,1], y=X_new[,2], fhat = fhat, sd = sqrt(var))
# ggplot(df, aes(x, y, fill= fhat)) +
#   geom_tile()

png("figures/pressure.png",width=800,height=600)
# png("figures/temp.png",width=800,height=600)
par(mfrow = c(1, 2))
cols <- heat.colors(100)
image(x_lat, x_long, matrix(fhat, ncol = n_grid), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 2], X[, 1])
image(x_lat, x_long, matrix(sqrt(var), ncol = n_grid), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 2], X[, 1])
dev.off()

