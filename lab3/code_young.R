# Question1

#Importing data
data1 = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab3/rainfall.dat", header = FALSE)

#1-1. Normal Model

## defining function for simulation of mu and sigmasq
library(extraDistr)
joint_simulation = function(data, niteration){
  
  
  ### initial nu is set to be 1. tau value is from the data.
  ### when scale parameter is fixed, the tail of the distribution is heavier
  ### as the degrees of freedom(nu) is smaller. Since the true sigmasq is 
  ### unknown, distribution with heavy tail is used to represent its uncertainty
  sigmasq = rinvchisq(1, nu=1, tau= var(data[,1]))
  init_mu = mean(data[,1])
  init_sigmasq = var(data[,1])
  
  n = dim(data)[1]

  result_store = data.frame(mu=init_mu, sigmasq = init_sigmasq)
  
    
  for (i in 1:niteration){
    w = (n/sigmasq) / ( (n/sigmasq) + (1/init_sigmasq) )
    tausq = 1 / ( (n/sigmasq) + (1/init_sigmasq) )
    new_mu = rnorm(1, mean = (w*(mean(data[,1]))) + ((1-w)*init_mu), sd = sqrt(tausq) )
    scaling = sum((data[,1] - new_mu)^2) / (n-1)
    new_sigmasq = rinvchisq(1, nu = n-1, tau = scaling)
    new = data.frame(mu = new_mu, sigmasq = new_sigmasq)
    result_store = rbind(result_store, new)
    init_mu = new_mu
    init_sigmasq = result_store$sigmasq[i]
    
  }
  
  return(result_store)
}

### simulation and plotting
set.seed(12345)
simulation_result = joint_simulation(data = data1, niteration = 1000)
### raw data plotting
plot(simulation_result$mu)
hist(simulation_result$mu)
plot(simulation_result$sigmasq)
hist(simulation_result$sigmasq, breaks = 100)
#### simulated mu values show somewhat similar distribution to normal distribution,
#### although it is a little bit left skewed.

### cumulative estimate plotting
cum_mu = cumsum(simulation_result$mu)
for (i in 1:length(cum_mu)){
  cum_mu[i] = cum_mu[i]/i
}
plot(cum_mu, main= "Trace plot of MCMC iteration of mu", xlab="MCMC iteration", ylab="Cummulative Estimate")

cum_sigmasq = cumsum(simulation_result$sigmasq)
for (i in 1:length(cum_sigmasq)) {
  cum_sigmasq[i] = cum_sigmasq[i]/i
}
plot(cum_sigmasq, main= "Trace plot of MCMC iteration of sigma^2", xlab="MCMC iteration", ylab="Cummulative Estimate")

### it is possible to see the cummulative estimate of the parameters are converging to certain values after burn-in period.



# 1-2. Mixture Normal

hist(data1[,1], breaks = 50, main = "Histogram of the Data", xlab = "Precipitation")
## assume that the data follows two-component mixture of normals model: one with
## low mean, and the other one with high mean

## change the form of data into matrix for more convenient use
#data(faithful)
#rawData <- faithful
x <- as.matrix(data1[,1])

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- c(20, 160) # Prior mean of mu by seeing the histogram of the data
tau2Prior <- rep(10,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(3,nComp) # degrees of freedom for prior on sigma2, set it as small value to make prior distribution of sigmasq vague

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green")
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))

storage_mu = data.frame(mu1 = as.numeric(), mu2 = as.numeric(),stringsAsFactors = FALSE)
storage_sigmasq = data.frame( sigmasq1 = as.numeric(), sigmasq2 = as.numeric(), stringsAsFactors = FALSE)

for (k in 1:nIter){
  #message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  #print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  new_mu = data.frame(mu1 = mu[1], mu2 = mu[2])
  storage_mu = rbind(storage_mu, new_mu)
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  new_sigma2 = data.frame(sigmasq1 = sigma2[1], sigmasq2 = sigma2[2])
  storage_sigmasq = rbind(storage_sigmasq, new_sigma2)
  
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    #hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      #lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    #lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    #legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
     #      col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted mixture density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 1, col = "red")
#lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density"), col=c("black","red"), lwd = 2)

#########################    Helper functions    ##############################################

### traceplot of mu1 and mu2
cum_mu1 = cumsum(storage_mu$mu1)
cum_mu2 = cumsum(storage_mu$mu2)
cum_sigmasq1 = cumsum(storage_sigmasq$sigmasq1)
cum_sigmasq2 = cumsum(storage_sigmasq$sigmasq2)
for (i in 1:nIter){
  cum_mu1[i] = cum_mu1[i]/i
  cum_mu2[i] = cum_mu2[i]/i
  cum_sigmasq1[i] = cum_sigmasq1[i]/i
  cum_sigmasq2[i] = cum_sigmasq2[i]/i
}

par(mfrow = c(1,2))
plot(cum_mu1, type = "l", main="Traceplot of mu1", xlab="MCMC Iteration")
plot(cum_mu2, type = "l", main ="Traceplot of mu2", xlab = "MCMC Iteration")
plot(cum_sigmasq1, type = "l", main="Traceplot of sigmasq1", xlab="MCMC Iteration")
plot(cum_sigmasq2, type = "l", main ="Traceplot of sigmasq2", xlab = "MCMC Iteration")



## 1-3.
### histogram of the given data
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Model Comparison")
### adding
#### density of mixture model using the mean of simulated mus and sigmasqs
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 1, col = "red")
#### density of normal model using the mean of simulated mu and sigmasq
lines(xGrid, dnorm(xGrid, mean=cum_mu[length(cum_mu)], sd = sqrt(cum_sigmasq[length(cum_sigmasq)])), type="l", lwd=2, lty=1, col="blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)



# Question2
## importing data & set variables
#data2 = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab3/eBayNumberOfBidderData.dat", header=TRUE)
data2 = read.table("lab3/eBayNumberOfBidderData.dat", header=TRUE)
Y = as.matrix(data2$nBids)
X = as.matrix(data2[,-1])

##2-1.
poisson_model = glm(nBids~0+., data = data2, family = "poisson")
summary(poisson_model)
poisson_model$coefficients
### significant covariates are:
#### Constant
#### VerifyID
#### Sealed
#### MajBlem
#### LogBook
#### MinBidShare ------------ these covariates have p-value less than 0.02

##2-2.
prior_cov = cov(X)


### basic setups
library(mvtnorm)
covariates = names(data2[,-1])
n_param = ncol(X)

### hyperparameters for priors
mu0 = rep(0, times = n_param)
sigma0 = 100*solve(t(X)%*%X)

### defining log_posterior functions used to find the mode of the distributions
minuslog_post = function(betas, x, y) {
  
  log_likelihood=rep(0,length(y))
  for(i in 1:length(y)){
    log_likelihood[i]=y[i]*betas%*%x[i,]-exp(betas%*%x[i,])-log(factorial(y[i]))
  }
  
  
  logprior = dmvnorm(betas, mean = as.matrix(mu0), sigma0, log=TRUE)
  # since it is log posterior, it is proportional to sum of loglikli&logprior
  # since optim function looks for minimum, minus posterior dist should be used
  return(-sum(log_likelihood)-logprior)
}

### objective is to find posterior distribution of betas.
### since the objective is to find the maximum, set fn as minus posterior
initials = rep(0, n_param)
optimum_result = optim(initials, minuslog_post, X, Y, gr = NULL, method = "BFGS", hessian = TRUE)

post_betas_mode = optimum_result$par  ### these turn out to be 0s... 
post_cov = -solve(optimum_result$hessian)
names(post_betas_mode) = covariates
colnames(post_cov) = covariates
rownames(post_cov) = covariates


## 2-3
library(mvtnorm)
RMW_sample = function(prev_betas, tuning, sigma, logposterior, ...) {
  proposal = rmvnorm(1, mean = prev_betas, sigma = tuning*sigma)
  alpha = min(1, logposterior(proposal, x, y)/logposterior(prev_betas, x, y))
  u = runif(1)
  if (u < alpha) {return(proposal)}
  else {return(prev_betas)}
}

example = matrix(NA, nrow = 10, ncol=n_param)
for (i in 1:10){
  example[i,] = RMW_sample(post_betas_mode, tuning = 1, sigma = post_cov, logposterior = minuslog_post, X, Y)
}
