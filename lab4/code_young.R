#1-a
ar1_sim = function(mu, phi, sigmasq, times) {
  x = mu
  res = c(x)
  for (i in 1:(times-1)){
    eps = rnorm(1, mean = 0, sd = sqrt(sigmasq))
    x = mu+(phi*(x-mu)) + eps
    res = append(res, x)
  }
  return(res)
}


mu = 10
sigmasq = 2
times = 200

phi_seq = c(-0.5, 0, 0.25,1)
res_mat = matrix(NA, ncol = times, nrow = length(phi_seq))

for (phi in phi_seq){
  res = ar1_sim(mu, phi, sigmasq, times)
  res_mat[which(phi_seq == phi),] = res
}

plot(res_mat[1,], type="l", col=1, main = "AR1 Plot")

for (i in 2:length(phi_seq)){
  lines(res_mat[i, ], col = i)
}

legend("bottomleft", legend = c("phi = -0.5", "phi = 0", "phi = 0.25", "phi=1"), 
       col = 1:4, lty = rep(1, times = 4))

## it seems as the absolute value of phi increases, oscillation range also increases


#1-b
library(rstan)
data_3 = ar1_sim(mu, phi=0.3, sigmasq, times)
data_95 = ar1_sim(mu, phi= 0.95, sigmasq, times)
sigma = sqrt(sigmasq)
N = length(data_3)

## constructin AR1 model for stan
StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigmasq;
}
model {
  for (n in 2:N)
    y[n] ~ normal(alpha + beta * (y[n-1]-alpha), sqrt(sigmasq));
}
'

data3 = list(N=N, y=data_3)
burnin=500
niter = 2000
fit3 = stan(model_code = StanModel, data= data3, warmup = burnin,
           iter = niter, chains = 4)
print(fit3, digits_summary = 3)
## mean, 95% CI, n_eff here. possible to estimate the true values
postDraws3 = extract(fit3)

data95 = list(N=N, y = data_95)
fit95 = stan(model_code = StanModel, data= data95, warmup = burnin,
             iter = niter, chains = 4)
print(fit95, digits_summary = 3)
## mean, CI, n_eff here. possible to estimate the true values
postDraws95 = extract(fit95)

## convergence of the samplers
### 0.3  =============> all converging
cum_mu3 = cumsum(postDraws3$alpha)
for (i in 1:length(cum_mu3)){
  cum_mu3[i] = cum_mu3[i]/i
}
plot(cum_mu3)
cum_phi3 = cumsum(postDraws3$beta)
for (i in 1:length(cum_phi3)){
  cum_phi3[i] = cum_phi3[i]/i
}
plot(cum_phi3)
cum_sigmasq3 = cumsum(postDraws3$sigmasq)
for(i in 1:length(cum_sigmasq3)){
  cum_sigmasq3[i] = cum_sigmasq3[i]/i
}
plot(cum_sigmasq3)
pairs(fit3)

### 0.95 =====================> all converging
cum_mu95 = cumsum(postDraws95$alpha)
for (i in 1:length(cum_mu95)){
  cum_mu95[i] = cum_mu95[i]/i
}
plot(cum_mu95)
cum_phi95 = cumsum(postDraws95$beta)
for (i in 1:length(cum_phi95)){
  cum_phi95[i] = cum_phi95[i]/i
}
plot(cum_phi95)
cum_sigmasq95 = cumsum(postDraws95$sigmasq)
for(i in 1:length(cum_sigmasq95)){
  cum_sigmasq95[i] = cum_sigmasq95[i]/i
}
plot(cum_sigmasq95)
pairs(fit95)
#### the joint posterior of mu and phi are concentrated around parameters' mean.
#### the joint posterior distribution seems normally distributed
#### since we know the data is simulated using normal distribution, it is possible to conclude that
#### the model constructed for sampling described the data well.



# 1-c.
#data = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab4/campy.dat", header = TRUE)
## constructin model for stan
"Pois_StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real<lower=0> mu; //because the parameter of poisson dist. must be positive
  real<lower = -1, upper = 1> phi;
  real<lower=0> sigmasq;
  vector[N] x;
}
model {
  mu~ normal(10,2);
  phi~ normal(0,100);
  sigmasq ~ scaled_inv_chi_square(1,2);
  x[1]~normal(mu, sqrt(sigmasq));
  x[2:N]~normal(mu + phi * (x[1:N-1]-mu), sqrt(sigmasq));
  y[1:N]~poisson(exp(x[1:N]));
}
'

input_data = list(N=dim(data)[1], y=data$c)
burnin=500
niter = 2000
fit = stan(model_code = Pois_StanModel, data= input_data, warmup = burnin,
            iter = niter, chains = 4)
"

data = read.table("~/Bayesian Learning/Bayesian-Learning-Lab/lab4/campy.dat", header = TRUE)
data = list(N=length(data$c), ct=data$c)
StanModel ='
data {
  int<lower=0> N;
  int ct[N];
}
parameters {
  real mu;
  real <lower=-1,upper=1>phi;
  real<lower=0> sigma;
  real xt[N];
}
model {
  mu~ normal(1,100);
  phi~ normal(1,100);
  sigma ~ scaled_inv_chi_square(1,2);
  xt[1]~normal(mu,1);
  ct[1]~poisson(exp(mu));
  for (n in 2:N){
    xt[n]~normal(mu + phi * (xt[n-1]-mu), sqrt(sigma));
    ct[n] ~ poisson(exp(xt[n]));
  }
}
'

burnin = 1000
niter = 2000
fit1 = stan(model_code=StanModel,data=data, warmup=burnin,iter=niter,chains=4)
print(fit1, digits_summary = 3)
postDraws = extract(fit1)
post_xt = postDraws$xt

means = c()
lowers = c()
uppers = c()
n = 4000 #number of samples

for (i in 1:ncol(post_xt)){
  means = append(means, mean(post_xt[,i]))
  temp = sort(post_xt[,i])
  lowers = append(lowers, temp[n*0.025])
  uppers = append(uppers, temp[n*0.975])
}

plot(data$ct,type="l", main = "Comparing Data and Simulation")
lines(exp(means), col = "red")
lines(exp(lowers), col = "blue")
lines(exp(uppers), col = "blue")
legend("topright", legend = c("Data", "Theta-mean", "95% C.I,"),
       col = c("black", "red", "blue"), lty = c(1,1,1))
