# beta-binomial model of allele frequencies
p <- seq(0, 1, 0.01)
k <- 40
n <- 200
alpha <- k+1
beta <- n-k+1

y <- dbeta(p, shape1=alpha, shape2=beta)

plot(x=p, y=y, ylab="Posterior density" , xlab="Population frequency of T", type="l")

qbeta(p=c(0.025,0.25,0.5,0.75,0.975), shape1=alpha, shape2=beta)

# with only 10 individuals
k <- 4
n <- 20

alpha <- k+1
beta <- n-k+1

y <- dbeta(p, shape1=alpha, shape2=beta)

plot(x=p, y=y, ylab="Posterior density" , xlab="Population frequency of T", type="l")

qbeta(p=c(0.025,0.25,0.5,0.75,0.975), shape1=alpha, shape2=beta)
# with informative prior

p <- seq(0, 1, 0.01)

k <- 40
n <- 200
alpha <- k+0.1
beta <- n-k+0.1
y <- dbeta(p, shape1=alpha, shape2=beta)
plot(x=p, y=y, ylab="Posterior density" , xlab="Population frequency of T", type="l")
qbeta(p=c(0.025,0.25,0.5,0.75,0.975), shape1=alpha, shape2=beta)

k <- 4
n <- 20
alpha <- k+0.1
beta <- n-k+0.1
y <- dbeta(p, shape1=alpha, shape2=beta)
points(x=p, y=y, type="l", lty=2)
qbeta(p=c(0.025,0.25,0.5,0.75,0.975), shape1=alpha, shape2=beta)

# bayes factors

th <- 0.1

p <- seq(0, 1, 0.01)
k <- 40
n <- 200

percs <- matrix(NA, nrow=3, ncol=5)
alpha <- k+1
beta <- n-k+1
y <- dbeta(p, shape1=alpha, shape2=beta)
plot(x=p, y=y, ylab="Posterior density" , xlab="Population frequency of T", type="l")
percs[1,1:4] <- c(qbeta(p=c(0.025,0.5,0.975), shape1=alpha, shape2=beta), pbeta(th, shape1=alpha, shape2=beta, lower.tail=FALSE) )

post1 <- pbeta(th, shape1=alpha, shape2=beta, lower.tail=FALSE)
prior1 <- pbeta(th, shape1=1, shape2=1, lower.tail=FALSE)
percs[1,5] <- (post1-(1-post1))/(prior1/(1-prior1))

alpha <- k+0.1
beta <- n-k+0.1
y <- dbeta(p, shape1=alpha, shape2=beta)
points(x=p, y=y, type="l", lty=2)
percs[2,1:4] <- c(qbeta(p=c(0.025,0.5,0.975), shape1=alpha, shape2=beta), pbeta(th, shape1=alpha, shape2=beta, lower.tail=FALSE) )


post1 <- pbeta(th, shape1=alpha, shape2=beta, lower.tail=FALSE)
prior1 <- pbeta(th, shape1=0.1, shape2=0.1, lower.tail=FALSE)
percs[2,5] <- (post1-(1-post1))/(prior1/(1-prior1))

alpha <- k+2
beta <- n-k+2
y <- dbeta(p, shape1=alpha, shape2=beta)
points(x=p, y=y, type="l", lty=3)
percs[3,1:4] <- c(qbeta(p=c(0.025,0.5,0.975), shape1=alpha, shape2=beta), pbeta(0.1, shape1=alpha, shape2=beta, lower.tail=FALSE) )

post1 <- pbeta(th, shape1=alpha, shape2=beta, lower.tail=FALSE)
prior1 <- pbeta(th, shape1=2, shape2=2, lower.tail=FALSE)
percs[3,5] <- (post1-(1-post1))/(prior1/(1-prior1))

