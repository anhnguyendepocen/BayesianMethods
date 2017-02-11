## Approximate Bayesian Computation (ABC)

# Rejection algorithm
N <- 1e5
y <- 4
simulate <- function(param) rpois(n=1, lambda=param)

thetas <- c()
while (length(thetas) <= N) {

	# 1. draw from prior (discrete, bounded, uniform)
	theta <- sample(0:10, 1)

	# 2. simulate observations
	ysim <- simulate(theta)

	# 3. accept/reject
	if (ysim == y) thetas <- c(thetas, theta)

}
hist(thetas)
quantile(thetas, c(0.025,0.25,0.5,0.75,0.975))

# Rejection algorithm in the continuous case
N <- 1e4
y <- 91.3514
epsilon <- 1e-1

# again this is supposed to be an unknown function
simulate <- function(param) rnorm(n=1, mean=param, sd=sqrt(10))

# euclidean distance
rho <- function(x,y) sqrt((x-y)^2)

thetas <- c()
while (length(thetas) <= N) {

	# 1. draw from prior (continuous, bounded, uniform)
	theta <- runif(1, min=50, max=150)

	# 2. simulate observations
	ysim <- simulate(theta)

	# 3. accept/reject
	if (rho(ysim,y)<=epsilon) thetas <- c(thetas, theta)

	cat(theta, ysim, rho(ysim,y), (rho(ysim,y)<=epsilon), length(thetas), "\n")
	Sys.sleep(2)
}
hist(thetas)
quantile(thetas, c(0.025,0.25,0.5,0.75,0.975))

# Rejection algorithm with proportions of simulations to accept
N <- 1e4
Y <- c(91.34, 89.21, 88.98)
th <- 0.05

simulate <- function(param) rnorm(n=1, mean=param, sd=sqrt(10))

# distance function
rho <- function(x,y) sum(sqrt((x-y)^2))/length(y)

thetas <- distances <- c()
for (i in 1:N) {

	# 1. draw from prior (continuous, bounded, uniform)
	theta <- 0
	while (theta<80 | theta>110) {
		theta <- rnorm(1, mean=90, sd=sqrt(20))
    }
	thetas <- c(thetas, theta)

	# 2. simulate observations
	ysim <- simulate(theta)

	# 3. calculate and retain distances
	distances <- c(distances, rho(ysim,Y))

}

# prior / parameter
hist(thetas)
quantile(thetas)

# relationship between parameter and distance
head(cbind(thetas,distances))

# distances
hist(distances)
accepted <- which(rank(distances, ties.method="random")/length(distances)<=th)
range(distances[accepted])
abline(v=max(distances[accepted]), lty=2)

# plot prior/posterior
par(mfrow=c(2,1))
hist(thetas, xlim=c(80,110), main="Prior")
hist(thetas[accepted], xlim=c(80,110), main="Posterior")
quantile(thetas[accepted])

# 95% HPD
library(coda)
HPDinterval(as.mcmc(thetas[accepted]), prob=0.95)
	
	

	
