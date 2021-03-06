# Question 1.
## 1-a.
s = 5
n = 20
f = n-s
alpha0 = 2
beta0 = 2
post_alpha = alpha0+s
post_beta = beta0+f
real_mean = post_alpha / (post_alpha+post_beta)
real_sd = sqrt( (post_alpha*post_beta) / (((post_alpha+post_beta)^2) *(post_alpha+post_beta+1)) )

means = c()
sds = c()
sequence = seq(from=10, to=10000, by=10)

set.seed(12345)
for (i in sequence){
  points = rbeta(n=i, shape1=post_alpha, shape2=post_beta)
  new_mean = mean(points)
  new_sd = sd(points)
  means = append(means, new_mean)
  sds = append(sds, new_sd)
}

plot(x = sequence, y = means, main = "Mean of the Posterior Distribution", xlab= "Number of Drawn Random Numbers", ylab = "Mean Value of the Drawn Random Numbers")
abline(h = real_mean, col="red")

plot(x = sequence, y = sds, main="Standard Deviation of the Posterior Distribution", xlab = "Number of Drqwn Random Numbers", ylab = "SD of the Drqwn Random Numbers")
abline(h = real_sd, col = "red")


## 1-b.
NDraw = 10000
set.seed(12345)
drawns = rbeta(n=NDraw, shape1 = post_alpha, shape2 = post_beta)
prob = sum(drawns > 0.3)/ NDraw
set.seed(12345)
true_prob3 = 1-pbeta(0.3, shape1 = post_alpha, shape2 = post_beta)
data.frame(Post_probability = prob, Exact_value = true_prob3, Difference_in_Percent = abs(prob - true_prob3)*100/true_prob3)
### difference in percentage is less than 0.5% ---- two values are very close to each other

## 1-c.
log_odds = log(drawns/(1-drawns))
hist(log_odds, xlab="log_odd values")


# Question 2

y = c(44,25,45,52, 30, 63, 19, 50, 34, 67)
log_y = log(y)
n = length(y)
mu = 3.7
tausq = sum((log_y - mu)^2)/n

## 2-a.

library(extraDistr)
set.seed(12345)
sigmas = rinvchisq(n=10000, nu=n, tau = tausq)
plot(x=1:10000, y=sigmas, main="Comparing Simulated Values and Real Value", xlab="Trial", ylab="Simulated Sigma-Squared")
theoretical_val=n*tausq/(n-2)
abline(h=theoretical_val, col="red")
data.frame(Simulation_Mean=mean(sigmas), Theoretical_Sigmasq = theoretical_val, Difference_in_Percentage = abs(mean(sigmas)-theoretical_val)*100/theoretical_val)
### difference between simulated mean and theoretical value in percentage is less than 1% ... very close values


## 2-b.
g_coeffs = (2*pnorm(sigmas/sqrt(2), mean=0, sd=1)) -1
hist(g_coeffs, main="Histogram of Gini Coefficients from Simulation", xlab="Simulated Gini-coefficients")

## 2-c.
## 90% equal tail credible interval
g_coeffs_sorted = sort(g_coeffs)
a = g_coeffs_sorted[10000*0.05]
b = g_coeffs_sorted[10000*0.95]
cat("90% equal tail credible interval is:", "\n")
c(a,b)
## kernel density estimate to compute 90% HPD interval
d_est <- density(g_coeffs_sorted)
sim_dens <- cumsum(d_est$y)/sum(d_est$y)
lower <- which(sim_dens>=0.05)[1]
upper <- which(sim_dens>=0.95)[1]
d_est$x[c(lower,upper)]
##compare two intervals-- HPD has wider interval


# Question 3
degrees = c(40, 303, 326, 285, 296, 314, 20, 308, 299, 296)
radians = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
n3 = length(radians)
known_mu = 2.39
## 3-a.
post_dist = function(kappa){
  res = exp((kappa*sum(cos(radians-known_mu))) - kappa) / (besselI(kappa, nu=0)^n3)
  return(res)
}
integrated = integrate(post_dist, lower=0, upper=20)
kappa_range = seq(0, 20, by=0.05)
plot(kappa_range, post_dist(kappa_range)/integrated$value, main="Posterior Distribution", xlab="Kappa", ylab="Density")
lines(kappa_range, post_dist(kappa_range)/integrated$value, col="red", lwd=2)
legend("topright", legend = "Posterior Distribution", col="red", lty=1)

## 3-b.
# definition of mode: the most frequent observation (to be more accurate in the report)
# therefore, the maximum value of density in posterior distribution will be equivelant to mode
cat("Approximate mode of K is: ", "\n")
kappa_range[which.max(post_dist(kappa_range))]
## division using constant is not needed because such division does not affect the result