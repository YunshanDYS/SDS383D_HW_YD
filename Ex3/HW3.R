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




