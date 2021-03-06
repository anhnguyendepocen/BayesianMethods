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


