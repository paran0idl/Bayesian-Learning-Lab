# q1
library(ggplot2)
set.seed("12345")
s=5
f=15
n=20
a0=2
beta0=2
true_a=a0+s
true_beta=beta0+f
means=c()
stds=c()
for(i in seq(10,10000,10)){
  pos=rbeta(i,true_a,true_beta)  
  means=c(means,mean(pos))
  stds=c(stds,sd(pos))
}

true_mean=true_a/(true_a+true_beta)
true_var=(true_a*true_beta)/((true_a+true_beta+1)*(true_a+true_beta)^2)

ggplot()+geom_point(aes(x=seq(10,10000,10),y=means))+geom_hline(yintercept=true_mean,col="red")
ggplot()+geom_point(aes(x=seq(10,10000,10),y=stds))+geom_hline(yintercept = sqrt(true_var),col="red")

nDraws=10000
res=rbeta(nDraws,true_a,true_beta)
prob=length(res[res>0.3])/nDraws
1-pbeta(0.3,true_a,true_beta)

res=rbeta(nDraws,true_a,true_beta)
log_odds=log(pos/(1-pos))
hist(log_odds)
plot(density(log_odds))


mu=3.7
income=c(44,25, 45, 52, 30, 63, 19, 50, 34,67)
tau_sqr=sum((log(income)-mu)^2)/length(income)
# sigma=seq(0,10,0.1)
# plot((1/sigma)^(length(income))*exp(sum((log(income)-mu)^2)*(-1/2*sigma^2)))


library(extraDistr)

nDraws=10000
x=rchisq(1,nDraws)
sigma=nDraws*tau_sqr/x
set.seed("12345")
inv=rinvchisq(10000,10,tau_sqr)

mean(inv)
ggplot()+geom_point(aes(x=1:10000,y=inv))+geom_hline(yintercept = sigma,col="red")
plot(x=seq(1:10000),y=inv)


# Gini
set.seed(12345)
Gini=2*pnorm(inv/sqrt(2))-1
hist(Gini)



#equal tail
pos_at_005=qinvchisq(0.05,10,tau_sqr)
pos_at_095=qinvchisq(0.95,10,tau_sqr)

2*pnorm(pos_at_005/sqrt(2))-1
2*pnorm(pos_at_095/sqrt(2))-1

# HPD
d=density(Gini)
which.max(d$y)
d$y[38]
d$x[38]
plot(d)
gini_sd=sd(Gini)
gini_mean=0.07014142
gini_mean-1.95*gini_sd
gini_mean+1.95*gini_sd

# q3
mu=2.39
data=c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)

k=seq(0,15,0.1)

dis=function(k){
  return(exp(k*sum(cos(data-mu))-k)*besselI(k,0)^(-length(data)))
}
inte=integrate(dis,lower = 0,upper = 15)


plot(dis(k)/inte$value)
res=res/sum(res)*15
plot(res)
plot(density(res))
p_k=(1/besselI(k,0))^length(data)
plot(p_k)


