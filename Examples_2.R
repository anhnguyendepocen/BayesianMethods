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

