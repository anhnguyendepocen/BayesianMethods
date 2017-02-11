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
