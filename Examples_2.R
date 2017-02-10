# Bayesian inference

## Gamma posterior asymmetric, one-tailed prior
alpha <- 0.5
beta <- 1
theta <- seq(0, 20, 0.1)
prior <- dgamma(x=theta, shape=alpha, scale=beta)
y <- 1
posterior <- dgamma(x=theta, shape=y+alpha, scale=1/(1+1/beta))
plot(x=theta, y=posterior, xlab=expression(theta), ylab="Density", type="l")
lines(theta, prior, lty=3)

## Interval estimation
theta <- seq(0, 10, 0.05)
alpha <- 2
beta <- 1
posterior <- dgamma(x=theta, shape=alpha, scale=beta)
plot(x=theta, y=posterior, xlab=expression(theta), ylab="Posterior density", type="l")

abline(h=0.1, lty=3)

library(coda)
x <- rgamma(n=1e5, shape=alpha, scale=beta)
hpd <- HPDinterval(as.mcmc(x), prob=0.87)
abline(v=hpd, lty=1)

a <- 1-0.87
eqi <- quantile(x, probs=c( a/2, 1-(a/2) ) )
abline(v=eqi, lty=2)

