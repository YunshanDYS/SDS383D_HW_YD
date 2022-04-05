library("reshape2") 

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



