# Question 1
## 1-a.
### import file
data1 = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab2/TempLinkoping.txt", header = TRUE)
### setting given hyperparameters
mu0 = c(-10,100,-100)
ohm0 = diag(0.01, ncol=3, nrow=3)
nu0 = 4
sigmasq0 = 1

### simulate from conjugate prior & see if reasonable
plot(data1, main = "Actual Data Points and Prior Regression Curves")

library(mvtnorm)
library(extraDistr)
set.seed(12345)

for (i in 1:20){
  prior_sigmasq = rinvchisq(n=1, nu=nu0, tau = sigmasq0)
  prior_beta = rmvnorm(n=1, mean=mu0, sigma = (prior_sigmasq*solve(ohm0)))
  prior_reg = prior_beta[1,1] + prior_beta[1,2]*data1$time + prior_beta[1,3]*data1$time^2 + rnorm(n=1, mean=0, sd = sqrt(sigmasq_prior))
  lines(x = data1$time, y=prior_reg, col=i+1)
}
## the range of regression line varies a lot, but it generally follows the pattern of the given data. Thus the given hyperparameters could be utilized to produce prior distribution.

## 1-b.
### simulate joint posterior
n = dim(data1)[1]
x = data.frame(intercept = rep(1, times = n), time = data1$time, timesq = data1$time^2)
X = as.matrix(x)
Y = as.matrix(data1$temp)
beta_hat = solve(t(X)%*%X) %*% t(X) %*% Y

mu_n = solve((t(X)%*%X) + ohm0) %*% ((t(X) %*% X %*% beta_hat) + (ohm0%*%mu0))
ohm_n = (t(X)%*%X) + ohm0
nu_n = nu0+n
sigmasq_n = nu0*sigmasq0 + ((t(Y)%*%Y) + (t(mu0)%*%ohm0%*%mu0) - (t(mu_n)%*%ohm_n%*%mu_n))
sigmasq_n = sigmasq_n / nu_n

post_betas = data.frame(beta0 = as.numeric(), beta1 = as.numeric(), beta2 = as.numeric(), stringsAsFactors = FALSE)
post_sigmasqs = data.frame(sigmasq = as.numeric(), stringsAsFactors = FALSE)

for (i in 1:10000){
  post_sigmasq = rinvchisq(n=1, nu=nu_n, tau = sigmasq_n)
  new = data.frame(sigmasq = post_sigmasq)
  post_sigmasqs = rbind(post_sigmasqs, new)
  post_beta = rmvnorm(n=1, mean=mu_n, sigma = (post_sigmasq*solve(ohm_n)))
  new = data.frame(beta0 = post_beta[1], beta1 = post_beta[2], beta2 = post_beta[3])
  post_betas = rbind(post_betas, new)
}

### plot marginal post.dist. of each parameter
hist(post_sigmasqs$sigmasq, main = "Histogram of Posterior Sigma-Squared", xlab="Simulated Posterior Value")
hist(post_betas$beta0, main = "Histogram of Posterior Beta0", xlab="Simulated Posterior Value")
hist(post_betas$beta1, main = "Histogram of Posterior Beta1", xlab="Simulated Posterior Value")
hist(post_betas$beta2, main = "Histogram of Posterior Beta2", xlab="Simulated Posterior Value")

### scatterplot and overlaying curves
plot(data1, main = "Data Points and Posterior Regression Curves")
#### posterior median
betas = as.matrix(post_betas)
sigmas = as.matrix(post_sigmasqs)
y_hat = matrix(0, nrow = dim(data1)[1], ncol=nrow(betas))
set.seed(12345)
for (i in 1:nrow(betas)) {
  error = as.matrix(rnorm(n=nrow(y_hat), mean=0, sd = sqrt(sigmas[i,1])))
  y_hat[,i] = X %*% betas[i,] + error
}     # row values will represent predicted value of y for given time, using nth simulated beta and sigma values
y_med = c()
for (i in 1:nrow(y_hat)){
  med_value = median(y_hat[i,])
  y_med = append(y_med, med_value)
}
lines(x = data1$time, y = y_med, col = "red")
#### posterior 95% credible interval
y_low = c()
y_up = c()
for (i in 1:nrow(y_hat)){
  values = sort(y_hat[i,])
  y_low = append(y_low, values[10000*0.025])
  y_up = append(y_up, values[10000*0.975])
}
lines(x = data1$time, y = y_low, col = "green")
lines(x = data1$time, y = y_up, col = "blue")
legend("bottom", legend = c("Posterior Median", "Upper Bound of Credible Interval", "Lower Bound of Credible Interval"),
       col = c("red", "blue", "green"), lty = c(1,1,1))

### the 95% equal tailed credible interval includes most of the data.
### this does make sense because the regression is based on joint posterior distribution, and such
### joint posterior distribution is derived by using the data.
### if the prior belief was very strong and such belief had completely different pattern than that of data points,
### then the result here could be different. However, regression using assumed prior distribution on the first step of this simulation
### had concave pattern - thus it did not greatly interfere the likelihood while producing
### posterior distribution of the parameters.

## 1-c.
x_tilde = c()
for (i in 1:ncol(y_hat)){
  x_tilde = append(x_tilde, data1$time[which.max(y_hat[,i])])
} imp


x_tilde = x_tilde*365
hist(x_tilde, main="Distribution of x_tilde")

## 1-d.
### one solution for mitigating potential problem is to set different hyperprior\
### for prior distribution of betas.
### set mu0 = c(0,0,0,0,0,0,0); set ohm0 to a matrix that its inverse matrix's diagonal 
### value would be very small so that covariance matrix of the prior distribution
### of betas would consist small diagonal values --- make prior have stronger 
### impact on posterior

### when prior distribution is certain, likelihood based on data has smaller
### impact on deriving posterior distribution. in above case, prior belief
### that betas will be 0s are very strong, it is unlikely the resulting
### posterior distribution of beta will consist high values for all elements.


# Question 2
## import data
data2 = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab2/WomenWork.dat", header=TRUE)

## 2-a.
