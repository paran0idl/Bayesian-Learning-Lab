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

library(extraDistr)
set.seed("12345")
nDraws=10000
x=rchisq(1,nDraws)
sigma=nDraws*tau_sqr/x
inv=rinvchisq(10000,10000,tau_sqr)

#ggplot()+geom_point(aes(x=1:10000,y=inv))+geom_hline(yintercept = sigma,col="red")
plot(x=seq(1:10000),y=inv)


# Gini
Gini=2*pnorm(inv/sqrt(2))-1




#equal tail
test=qinvchisq(0.05,10000,tau_sqr)
test1=qinvchisq(0.95,10000,tau_sqr)

2*pnorm(test/sqrt(2))-1
2*pnorm(test1/sqrt(2))-1

# HPD
d=density(Gini)
gini_sd=sd(Gini)
gini_mean=0.08444
gini_mean-1.95*gini_sd
gini_mean+1.95*gini_sd

# q3
data=c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)

k=rexp(10000,1)
k=k[order(k)]
p_k=(1/besselI(k,0))^length(data)
plot(p_k)


