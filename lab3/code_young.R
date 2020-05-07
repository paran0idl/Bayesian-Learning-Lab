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
data(faithful)
rawData <- faithful
x <- as.matrix(rawData['eruptions'])

# Model options
nComp <- 4    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(0,nComp) # Prior mean of mu
tau2Prior <- rep(10,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
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


for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
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
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
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
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

#########################    Helper functions    ##############################################
