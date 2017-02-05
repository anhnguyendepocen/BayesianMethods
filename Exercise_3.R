
# Approximate Bayesian Computation

## PREPARATION

# open R

# load all R functions and data we need
source("functions.R")
load("polar.brown.sfs.Rdata")

ls()

# the file "polar.brown.sfs" includes the joint (2 dimensions) site frequency spectrum (SFS) between polar bears (on the rows) and brown bears (on the columns)
# if you want to see this file type `cat polar.brown.sfs` in your terminal

# we can plot it
plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="2D-SFS")

# each population has 2n+1 entries in its spectrum, with n being the number of individuals
# the number of chromosomes can be retrieved as
nChroms.polar < nrow(polar.brown.sfs)-1
nChroms.polar
nChroms.brown <- ncol(polar.brown.sfs)-1
nChroms.brown

# the number of analysed sites (here all the sites polymorphic) is simply the sum of all entries in the SFS
nrSites <- sum(polar.brown.sfs, na.rm=T)
nrSites
# we will condition the simulations to generate this number of sites for each repetition

# we can calculate the observed summary statistics
obsSummaryStats <- calcSummaryStats(polar.brown.sfs)

# the parameters we want to estimate are the divergence time between polar and brown bears (T) and the migration rate (M)

# we will use the `abc` package and the `abc` function
library(abc)
?abc

# as you can see we need 3 objects:
# target: a vector of the observed summary statistics.
# param: a vector, matrix or data frame of the simulated parameter values.
# sumstat: a vector, matrix or data frame of the simulated summary statistics.

# we already have 'target' as it is the vector of observed summary statistics called 'obsSummaryStatistics'

## TASK 1

# our first aim is to performs N simulations of data by drawing from prior distributions of T and M and record (separately) the drawn values and the corresponding summary statistics

# we can define how many simulations we want to perform (ideally a lot)
nrSimul <- 1e4

# then we define the prior distribution of our parameters to be estimated, divergence time T and migration rate M
# you can use any distribution you find suitable (e.g. uniform, normal, beta); however you may want to consider that a reasonable range of values for T is between 200k and 800k years ago; the migration rate is scaled by the reference population size so you can consider a reasonble range of M being between 0 (included) and 5.

# the function to simulate data (specifically the site frequency spectrum) given values of T and and M is 'simulate'
simulate
# and it takes as parameters: T, M plus some fixed value of how many sites, the directory for ms program and the text file in output.
# as an example, assuming T=200k and M=1 the command to simulate data and calculate summary statistics is the following:
# first, set the directory wher you installed "ms" software
msDir <- "~/Software/msdir/ms" # this is my specific case, yours could be different
# second, set the name for the output text file
fout <- "ms.txt" # leave it like here
# then we can simulate data:
simulate(T=2e5, M=1, nrSites, msDir, fout)
# and finally calculate the summary statistics for this simulation
simulatedSFS <- fromMStoSFS(fout, nrSites, nChroms.polar, nChroms.brown)
calcSummaryStats(simulatedSFS)
# you can even plot it
plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="simulated 2D-SFS")

# Question: based on the observed summary statistics compared to these ones using {T=2e5, M=1}, can you make some considerations on the most likely values of T or M (higher or lower)?

# we are now ready to perform the simulations; for each simulation we need to retain both the drawn parameters and the corresponding summary statistics
params <- matrix(NA, nrow=nrSimul, ncol=2)
colnames(params) <- c("T","M")
stats <- matrix(NA, nrow=nrSimul, ncol=length(obsSummaryStats))

# iterate across all simulations
for (i in 1:nrSimul) 
{

	# print
	cat("\n",i, "/", nrSimul)

	# draw parameters from prior distributions
	# e.g if uniform
	T_random <- runif(1, min=2e5, max=8e5)
	M_random <- runif(1, min=0, max=5)
	params[i,] <- c(T_random, M_random)

	# simulate data with these values
	simulate(T=T_random, M=M_random, nrSites, msDir, fout)
	# calculate and record summary statistics
	stats[i,] <- calcSummaryStats(fromMStoSFS(fout, nrSites, nChroms.polar, nChroms.brown))

}

## TASK 2

# inspect simulations
# range?
# correlation between summary statistics?
# correlation between summary stats and params?

# which to retain? which are informative?

# scale them (with the obs stats) to have mean=0 and sd=1

## TASK 3

# inference, using rejection or loc-linear

# plot posteriors (and priors)

# calculate bayes factors for each bin?

## TASK 4

# test model with no migration vs with migration






















