# Examples chapter 2

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

