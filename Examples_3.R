# Bayesian computation

# beta-binomial model of allele frequencies
p <- seq(0, 1, 0.01)

# 1000 chromosomes (20 derived alleles)
k <- 200
n <- 1000
alpha <- k+1
beta <- n-k+1
plot(x=p, y=dbeta(p, shape1=alpha, shape2=beta), ylab="Posterior density" , xlab="Population frequency of T", type="l")

# 100 chromosomes
k <- 20
n <- 100
alpha <- k+1
beta <- n-k+1
points(x=p, y=dbeta(p, shape1=alpha, shape2=beta), type="l", lty=2)

# 10 chromosomes
k <- 2
n <- 10
alpha <- k+1
beta <- n-k+1
points(x=p, y=dbeta(p, shape1=alpha, shape2=beta), type="l", lty=3)

legend("topright", legend=c(1e3,1e2,1e1), lty=1:3)

# beta-binomial model of allele frequencies
p <- seq(0, 1, 0.01)

# 100 chromosomes (20 derived alleles)
k <- 20
n <- 100

# exact posterior with flat prior
alpha <- k+1
beta <- n-k+1
plot(x=p, y=dbeta(p, shape1=alpha, shape2=beta), ylab="Density" , xlab="Population frequency", type="l")

# normal approximation
thetaHat <- k/n
var <- thetaHat*(1-thetaHat)/n
points(x=p, y=dnorm(p, mean=thetaHat, sd=sqrt(var)), type="l", lty=2)

legend("topright", c("exact","approx"), lty=1:2)

# data
k <- 20
n <- 100

# exact posterior with flat prior
alpha <- k+1
beta <- n-k+1
exact <- rbeta(1e5, shape1=alpha, shape2=beta)
quantile(exact)

# normal approximation
thetaHat <- k/n
var <- thetaHat*(1-thetaHat)/n
approx <- rnorm(1e5, mean=thetaHat, sd=sqrt(var))
quantile(approx)

# calculate density from true posterior
calcTrueDensity <- function(x) dbeta(x, 3, 10)

# g is a uniform prior, M is the maximum density value of the posterior (if known)
x <- seq(0,1,0.01)
epsilon <- 1e-3
M <- max(calcTrueDensity(x)) + epsilon

thetas <- c()

# we want N samples
N <- 1e4
rawDensity <- rbeta(N, 3, 10)

while (length(thetas) < N) {

	theta_j <- runif(1, 0, 1)
	U <- runif(1, 0, 1)
    
    if (M*U < calcTrueDensity(theta_j)) thetas <- c(thetas, theta_j) 

}

# check qq-plot
qqplot(rawDensity, thetas)
abline(0,1)
