# Examples 1 - Bayesian (thinking &) concepts

## normal/normal model

# prior
mu <- 2
tau <- 1
x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6), 
type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

# likelihood
y <- 6 
sigma <- 1
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2)

# posterior
B <- sigma^2/(sigma^2+tau^2)
postMean <- B*mu + (1-B)*y
postVar <- B*tau^2
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)

# skewer prior
mu <- 2
tau <- 0.5
x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,1), 
type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

# likelihood (obviously it does not change)
y <- 6 
sigma <- 1
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2)

# posterior with skewer prior
B <- sigma^2/(sigma^2+tau^2)
postMean <- B*mu + (1-B)*y
postVar <- B*tau^2
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)

# MAP with skewer prior
map <-  x[which.max(dnorm(x=x, mean=postMean, sd=sqrt(postVar)))]
cat("MAP:", map, "(", map-3*sqrt(postVar),",",map+3*sqrt(postVar),")\n")

# wider prior
mu <- 2
tau <- 2
x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.5), 
type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

# likelihood (obviously it does not change)
y <- 6 
sigma <- 1
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2)

# posterior with wider prior
B <- sigma^2/(sigma^2+tau^2)
postMean <- B*mu + (1-B)*y
postVar <- B*tau^2
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)

# MAP with wider prior
map <-  x[which.max(dnorm(x=x, mean=postMean, sd=sqrt(postVar)))]
cat("MAP:", map, "(", map-3*sqrt(postVar),",",map+3*sqrt(postVar),")\n")

# more observations (obviously the prior does not change)
mu <- 2
tau <- 1
x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,1), 
type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

# likelihood with more observations
y <- c(6, 5, 5.5, 6.5, 6)
n <- length(y)
sigma <- 1
points(x=x, y=dnorm(x=x, mean=mean(y), sd=sigma), type="l", lty=2)

# posterior with more observations
postMean <- ( (sigma^2/n)*mu + tau^2*mean(y) ) / ( (sigma^2/n)*mu + tau^2 )
postVar <- ( (sigma^2/n)*tau^2 ) / ( (sigma^2/n) + tau^2 )
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)

# MAP with more observations
map <-  x[which.max(dnorm(x=x, mean=postMean, sd=sqrt(postVar)))]
cat("MAP:", map, "(", map-3*sqrt(postVar),",",map+3*sqrt(postVar),")\n")

## monte carlo sampling

par(mfrow=c(3,1))

# posterior
x <- seq(2,8,0.01)
postMean <- 4.43
postVar <- 0.16
plot(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=1, ylab="Density", xlab=expression(theta), main=expression(p(theta~"|"~y)), ylim=c(0,1.2), xlim=c(2,8))

# sampling
y_sampled_1 <- rnorm(n=50, mean=postMean, sd=sqrt(postVar))
hist(y_sampled_1, breaks=20, freq=F, lty=2, col="grey", ylim=c(0,1.2), xlim=c(2,8), sub="50 samples", main=expression(p(theta~"|"~y)), xlab=expression(theta))
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=1)

# more sampling
y_sampled_2 <- rnorm(n=1e6, mean=postMean, sd=sqrt(postVar))
hist(y_sampled_2, breaks=20, freq=F, lty=2, col="grey", ylim=c(0,1.2), xlim=c(2,8), sub="1e6 samples", main=expression(p(theta~"|"~y)), xlab=expression(theta))
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=1, ylab="Density", xlab=expression(theta), main=expression(p(theta~"|"~y)), sub="1e6 samples")

## elicited prior
mu <- 88
sigma2 <- 10
x <- seq(80,110,0.1)
plot(x=x, dnorm(x=x, mean=mu, sd=sqrt(sigma2)),  
type="l", lty=1, ylab="Density", xlab=expression(theta),
main=expression(pi(theta)))

## Poisson
theta <- 4
y <- seq(0, 20, 1)
plot(x=y, dpois(x=y, lambda=theta), type="p", lty=1, xlab=expression(y), main=expression(theta~"="~4), ylab=expression(p(Y~"="~y~"|"~theta)))

## Gamma
alpha <- c(0.5,2,10)
beta <- c(0.5,2)
thetas <- seq(0, 20, 0.1)
mycolors <- topo.colors(6)
plot(x=thetas, dgamma(x=thetas, shape=alpha[1], scale=beta[1]), type="l", lty=1, xlab=expression(theta), main="Gamma", ylab=expression(pi(theta)), ylim=c(0,1.0), col=mycolors[1], lwd=2)
index <- 0
for (i in alpha) {
	for (j in beta) {
    	index <- index+1
    	points(x=thetas, dgamma(x=thetas, shape=i, scale=j), col=mycolors[index], ty="l", lwd=2)
    }
}
names <- cbind(rep(alpha,each=2),rep(beta))
legend(x="topright", legend=apply(FUN=paste, MAR=1, X=names, sep="", collapse=", "), col=mycolors, lty=1, lwd=2)

## Gamma posterior
alpha <- 3
beta <- 1
theta <- seq(0, 20, 0.1)
prior <- dgamma(x=theta, shape=alpha, scale=beta)
y <- 4
posterior <- dgamma(x=theta, shape=y+alpha, scale=1/(1+1/beta))
plot(x=theta, y=posterior, xlab=expression(theta), ylab="Density", type="l")
lines(theta, prior, lty=3)
postdraw <- rgamma(n=1e5, shape=y+alpha, scale=1/(1+1/beta))
histdraw <- hist(postdraw, breaks=20, plot=F)
lines(histdraw, lty=3, col="grey", freq=F)

## Uniform
theta <- seq(0, 100, 0.1)
prior <- dunif(x=theta, min=0, max=100)
plot(x=theta, y=prior, xlab=expression(theta), ylab="Density", type="l")

