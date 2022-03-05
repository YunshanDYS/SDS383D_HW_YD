## Ex 4
## data
data <- read.csv("mathtest.csv",header=TRUE)

library(dplyr)
y_school <- data %>%
  group_by(school) %>%
  summarize(bar_y = mean(mathscore), n_data = n())
## plot
jpeg("figures/bary_n.jpeg",width=600,height=800)
plot(y_school$bar_y,y_school$n_data,xlab="bar_y",ylab="n")
dev.off()