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
    scaling = (sigmasq + sum((data[,1] - new_mu)^2) / (n+1))
    new_sigmasq = rinvchisq(1, nu = 1+i, tau = scaling)
    new = data.frame(mu = new_mu, sigmasq = new_sigmasq)
    result_store = rbind(result_store, new)
    init_mu = new_mu
    init_sigmasq = result_store$sigmasq[i]
    
  }
  
  return(result_store)
}

### simulation and plotting
simulation_result = joint_simulation(data = data1, niteration = 1000)
plot(simulation_result$mu)
hist(simulation_result$mu)
plot(simulation_result$sigmasq, ylim = c(0,3000))
hist(simulation_result$sigmasq, xlim = c(0,3000), breaks = 1000)
#### simulated mu values show somewhat similar distribution to normal distribution,
#### although it is a little bit left skewed. no clear convergence -------- check code needed i think
##### one possible reason: maybe it is because simulated mu value was already very similar
##### to the real mu so the traceplot does not show clear convergence. standard deviation
##### of the simulated mu values are 0.206, which is quite small

#### simulated sigmasq values converges to certain value as the simulation is
#### repeated.



# 1-2. Mixture Normal

hist(data1[,1], breaks = 50, main = "Histogram of the Data", xlab = "Precipitation")
## assume that the data follows two-component mixture of normals model: one with
## low mean, and the other one with high mean

## change the form of data into matrix for more convenient use
data = as.matrix(data1)
ncomp = 2 # because it is of two-component