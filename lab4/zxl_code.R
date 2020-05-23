library(ggplot2)
set.seed("12345")
ar1_process=function(t,mu,phi,sigmasq){
  res=vector(length=t)
  res[1]=mu
  for(i in 2:t){
    res[i]=mu+phi*(res[i-1]-mu)+rnorm(1,0,sigmasq)
  }
  return(res)
}
mu=10
sigmasq=2
t=200

tmp=ar1_process(t,mu,0.5,sigmasq)
plot(x=seq(1,t),y=tmp)
#ggplot()+geom_line(aes(x=seq(1,t),y=tmp))

x1t=ar1_process(200,10,0.3,2)
y1t=ar1_process(200,10,0.95,2)
mean(x1t)
mean(y1t)
ggplot()+geom_line(aes(x=seq(1,t),y=y1t))

library(rstan)

StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real <lower=-1,upper=1>phi;
  real<lower=0> sigma;
}
model {
  mu~ normal(1,100);
  phi~ normal(1,100);
  sigma ~ scaled_inv_chi_square(1,2);
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1]-mu), sigma);
}
'

data1 = list(N=length(x1t), y=x1t)
data2 = list(N=length(y1t), y=y1t)
burnin = 1000
niter = 2000
fit1 = stan(model_code=StanModel,data=data1, warmup=burnin,iter=niter,chains=4)
fit2 = stan(model_code=StanModel,data=data2, warmup=burnin,iter=niter,chains=4)
# Print the fitted model
print(fit1,digits_summary=3)
print(fit2,digits_summary=3)
phi=extract(fit1)$phi
hist(phi,breaks=100)
mu=extract(fit1)$mu

# Extract posterior samples
postDraws <- extract(fit1)

# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$alpha,type="l")
plot(postDraws$beta,type="l")
plot(postDraws$sigma,type="l")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)


data=read.table("lab4/campy.dat",header = TRUE)
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
    xt[n]~normal(mu + phi * (xt[n-1]-mu), sigma);
    ct[n] ~ poisson(exp(xt[n]));
  }
}
'

plot(data$ct,type="l")
burnin = 1000
niter = 2000
fit1 = stan(model_code=StanModel,data=data, warmup=burnin,iter=niter,chains=4)


print(fit1,digits_summary=3)
postDraws <- extract(fit1)
res=colMeans(postDraws$xt)
plot(exp(res),type="l")

#
theta_data=list(N=length(res), ct=res)

StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real <lower=-1,upper=1>phi;
  real<lower=0> sigma;
}
model {
  mu~ normal(1,100);
  phi~ normal(1,100);
  sigma ~ scaled_inv_chi_square(1,2);
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1]-mu), sigma);
    //temp ~ normal(mu + phi * (y[n-1]-mu), sigma);
    
}
'

theta_data = list(N=length(res), y=res)
burnin = 1000
niter = 2000
theta_fit = stan(model_code=StanModel,data=theta_data, warmup=burnin,iter=niter,chains=4)
print(theta_fit,digits_summary=3)

test=ar1_process(140,3.484,0.926,0.173)

test=exp(test)
plot(rpois(140,test),type="l")