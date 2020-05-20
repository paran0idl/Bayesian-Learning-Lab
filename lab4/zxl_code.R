library(ggplot2)
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
#plot(x=seq(1,t),y=tmp)
ggplot()+geom_line(aes(x=seq(1,t),y=tmp))

x1t=ar1_process(200,10,0.3,2)
y1t=ar1_process(200,10,0.95,2)
ggplot()+geom_line(aes(x=seq(1,t),y=y1t))

library(rstan)
y=x1t
N=length(y)
mean(y1t)
var(y1t)
StanModel = '
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha~ normal(10,1); // Normal with mean 0, st.dev. 100
  beta~ normal(0,1); // Normal with mean 0, st.dev. 100
  sigma ~ scaled_inv_chi_square(1,5); // Scaled-inv-chi2 with nu 1, sigma 2
  for (n in 2:N)
    y[n] ~ normal(alpha + beta * (y[n-1]-alpha), sigma);
}
'




# library(rstan)
# y = c(4,5,6,4,0,2,5,3,8,6,10,8)
# N = length(y)


data = list(N=N, y=y)
burnin = 1000
niter = 2000
fit = stan(model_code=StanModel,data=data, warmup=burnin,iter=niter,chains=2)
# Print the fitted model
print(fit,digits_summary=3)
# Extract posterior samples
postDraws <- extract(fit)

# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws$alpha,type="l")
plot(postDraws$beta,type="l")
plot(postDraws$sigma,type="l")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)