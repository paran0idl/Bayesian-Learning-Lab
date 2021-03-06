set.seed("12345")
library(mvtnorm)
data=read.table("lab2/TempLinkoping.txt",header = TRUE)
x=data.frame(x1=rep(1,length(data$time)),x2=data$time,x3=data$time^2)
x=as.matrix(x)
y=data$temp
mu_0=c(-10,100,-100)
Omega_0=0.01*diag(1,3,3)
v_0=4
sigma_sq_0=1

# prior
library(extraDistr)
plot(data$temp)
for(i in 1:10){
  sigma_prior=rinvchisq(n=1, nu=v_0, tau = sigma_sq_0)
  beta_prior=rmvnorm(1,mean = mu_0,sigma = sigma_prior*solve(Omega_0))
  res=beta_prior[1,1]+beta_prior[1,2]*data$time+beta_prior[1,3]*data$time^2+rnorm(1,mean=0,sd=sqrt(sigma_prior))
  lines(res,col="red")
}
# all curves looks reasonable

# marginal posteriors

beta_hat=solve(t(x)%*%x)%*%t(x)%*%as.matrix(y)

mu_n=solve((t(x)%*%x+Omega_0))%*%(t(x)%*%x%*%beta_hat+Omega_0%*%mu_0)
Omega_n=t(x)%*%x+Omega_0
v_n=v_0+length(data$temp)
sigma_sq_n=v_0*sigma_sq_0+(t(y)%*%y+t(mu_0)%*%Omega_0%*%mu_0-t(mu_n)%*%Omega_n%*%mu_n)
sigma_sq_n=sigma_sq_n/v_n
sigma_sq=rinvchisq(1,v_n,sigma_sq_n)

betas=rmvnorm(10000,mean=mu_n,sigma=sigma_sq*solve(Omega_n))
epsilon=matrix(rnorm(10000,mean=0,sd=sqrt(sigma_sq)))
parms=data.frame(betas,epsilon)
# hist(betas[,1],breaks=100)
# hist(betas[,2],breaks=100)
# hist(betas[,3],breaks=100)

fn_time=function(parms,time,sigma_value){
  df=data.frame()
  for(i in 1:length(time)){
    res=parms[,1]+parms[,2]*time[i]+parms[,3]*time[i]^2+parms[,4]
    median_temp=median(res)
    res=sort(res)
    lower_va=res[length(res)*0.025]
    upper_va=res[length(res)*0.975]
    df=rbind(df,data.frame(median=median_temp,lower=lower_va,upper=upper_va))
  }
  return(df)
}
res=fn_time(parms,data$time,sigma_sq)

library(ggplot2)
ggplot()+geom_point(aes(x=data$time,y=data$temp))+
  geom_line(aes(x=data$time,y=res$median,colour="Median"))+
  geom_line(aes(x=data$time,y=res$lower,colour="Lower"))+
  geom_line(aes(x=data$time,y=res$upper,colour="Upper"))
  
# contains most of the points


# c
plot(density(parms$X2/(-2*parms$X3)))
# around 0.545

# d
# check slides




# q2
library(mvtnorm)
data=read.table("lab2/WomenWork.dat",header = TRUE)

# likelihood for Probit model 
# https://en.wikipedia.org/wiki/Probit_model

y=data$Work
x=data[,-1]
x=as.matrix(x)
co_count=dim(x)[2]

# log_likelihood_probit=function(y,x,betas){
#   return(sum(y*log(pnorm(x%*%betas))+(1-y)*log(1-pnorm(x%*%betas))))
# }

log_likehood_logistic=function(y,x,betas){
  numerator=exp(x%*%betas)
  return(log(prod(numerator^y/(1+numerator))))
}


log_prior_likelihood=function(tau,betas){
  sigmas=diag(10,co_count,co_count)
  mu=rep(0,co_count)
  log_likelihood=(log(det(sigmas))+t(betas-mu)%*%solve(sigmas)%*%(betas-mu)+co_count*log(2*pi))*0.5
  return(log_likelihood)
}

# for log posterior= log prior +log likelihood
log_posterior=function(betas,y,x,tau){
  #return(log_likelihood_probit(y,x,betas)+log_prior_likelihood(tau,betas))
  return(log_likehood_logistic(y,x,betas)+log_prior_likelihood(tau,betas))
}

starting_value=rep(0,co_count)
tau=10
op=optim(starting_value,log_posterior,gr=NULL,y,x,tau,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
op$par
glm_model=glm(Work ~ 0 + ., data = data, family = binomial)
glm_model$coefficients

posterior_mode=op$par
posterior_cov=-solve(op$hessian)
# res=x%*%posterior_mode
# res=as.numeric(res>0)
# 
# cft=table(res,data$Work)
# (cft[1,1]+cft[2,2])/sum(cft)


new_woman=c(1,10,8,10,1,40,1,1)
posterior_draws=rmvnorm(1000,mean=posterior_mode,sigma = posterior_cov)
new_woman_res=posterior_draws%*%new_woman

plot(density(new_woman_res))

p_work=c()
for(i in 1:10){
  posterior_draws=rmvnorm(1000,mean=posterior_mode,sigma = posterior_cov)
  new_woman_res=posterior_draws%*%new_woman
  new_woman_res=as.numeric(new_woman_res>0)
  work=length(new_woman_res[new_woman_res==1])
  p_work=c(p_work,work/1000)
}
p_work
# for binomial, how we get the probability for success ? averging or take one value ?
plot(density(rbinom(10000,10,0.3)))
?rbinom