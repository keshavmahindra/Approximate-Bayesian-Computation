rm(list=ls()) #First clear the R environment

observedData =  rbeta(50, 6, 3) #Generate the observed sample
observedSummary = c(mean(observedData), var(observedData)) #Calculate the sufficient (or summary) statistics.
 
# The summary here is used because simulation-based methods are typically more efficient 
# when the dimensionality of the data is low.
# In general, one has to check whether information is lost by such a reduction (sufficiency). 
# Here the mean and variance between them uniquely determine a and b, so they are sufficient statistics.
 
# Defining a stochastic model with beta output
# For convenience we are calculating the summary statistics in the same step
 
model <- function(par){
  simulatedData <- rbeta(50, par[1,1], par[1,2])
  simulatedSummary <- c(mean(simulatedData), var(simulatedData))
  return(simulatedSummary)
}

n = 10000 #sample size for the simulation - best not to set this very much larger if you use this code (I ran a 100,000 point simulation and it took a whole day to complete).
fit = data.frame(alpha = runif(n, 0,25), beta = runif(n, 0,10), 
summary1 = rep(NA, n), summary2 = rep(NA, n), distance = rep(NA, n))
 
#Using a for loop is super slow here, but good from the point of view of understanding the code. May want to vectorize this part.
for (i in 1:n){
  prediction <- model(fit[i,1:2])
  deviation = sqrt(sum((prediction- observedSummary)^2))
#A potential problem here is if the scale of one of the parameters is much higher than the other (for Beta(3,6) probably okay).
#If one parameter takes values much larger than the other then it will tend to dominate the distance - you end up fitting just that parameter component well in practice.
#Can be corrected for by normalization - i.e. standardize both values to have variance 1 *just for the purposes of the distance calculation*.
  fit[i,3:5] = c(prediction, deviation) #Store the simulated values *and their deviations*. NB we do not do *any* rejection yet.
}

plot(fit[fit[,5] < 0.5, 1:2], xlim = c(0,25), ylim = c(0,10), col = "lightgrey", main = "Accepted parameters for different values of epsilon")
points(fit[fit[,5] < 0.1, 1:2],  pch = 18, col = "gray")
points(fit[fit[,5] < 0.05, 1:2],  pch = 8, col = "red")
points(fit[fit[,5] < 0.025, 1:2],  pch = 20, col = "black")
points(fit[fit[,5] < 0.01, 1:2],  pch = 24, col = "white")
 
legend("bottomright", c("< 0.5", "< 0.1", "< 0.05", "< 0.025", "< 0.01"), pch = c(1,18,8,20,24), col = c("lightgrey", "gray", "red", "black", "white"))
 
#To see what the sampling looks like, we plot the target params as crosshairs!
abline(v = 6)
abline(h = 3)  
 
#Let's look at some summary stats and histograms of values for one choice of the rejection radius epsilon:
outcomes = fit[fit[,5] < 0.01, 1:2]
dim(outcomes)
summary(outcomes)
hist(outcomes[,1], freq = FALSE, main = 'Histogram of X', xlab = 'X')
hist(outcomes[,2], freq = FALSE, main = 'Histogram of X', xlab = 'X')

plot.new()
x = seq(0,1,0.00001)
mode = (6-1)/(6+3-2) #Analytic form of mode of a Beta distribution is (a-1)/(a+b-2)
modeHeight = dbeta(mode,6,3) #Work out height at the mode.
modeHeight
plot.window(xlim = c(0,1), ylim = c(0,(modeHeight*1.5)))
lines(x, dbeta(x,6,3), col="red")
lines(x, dbeta(x,mean(outcomes[,1]),mean(outcomes[,2])), col="blue")
outcomes2 = fit[fit[,5] < 0.02, 1:2]
dim(outcomes2)
summary(outcomes2)
lines(x, dbeta(x,mean(outcomes2[,1]),mean(outcomes2[,2])), col="green")
outcomes3 = fit[fit[,5] < 0.005, 1:2]
dim(outcomes3)
summary(outcomes3)
lines(x, dbeta(x,mean(outcomes3[,1]),mean(outcomes3[,2])), col="black")
outcomes4 = fit[fit[,5] < 0.0001, 1:2]
lines(x, dbeta(x,mean(outcomes4[,1]),mean(outcomes4[,2])), col="lightblue")
dim(outcomes4)
summary(outcomes4)
outcomes5 = fit[fit[,5] < 0.015, 1:2]
dim(outcomes5)
summary(outcomes5)
lines(x, dbeta(x,mean(outcomes5[,1]),mean(outcomes5[,2])), col="lightblue")
outcomes6 = fit[fit[,5] < 0.0075, 1:2]
dim(outcomes6)
summary(outcomes6)
lines(x, dbeta(x,mean(outcomes6[,1]),mean(outcomes6[,2])), col="darkcyan")
hist(outcomes6[,1], freq = FALSE, main = 'Histogram of X', xlab = 'X')
hist(outcomes6[,2], freq = FALSE, main = 'Histogram of X', xlab = 'X')
#Sweet spot seems to be around epsilon = 0.0075 - we can make a more fine-grained search in this region if necessary.
#Alternative ways of choosing epsilon (automatically) include cross-validation 
#(hold out some of the posterior sample in constructing the ABC model, choose the epsilon that gives smallest MSE on 
# the held out data.) and linear and non-linear regression techniques.
# You can use library(abc) to calculate the acceptance radius epsilon in these ways. 