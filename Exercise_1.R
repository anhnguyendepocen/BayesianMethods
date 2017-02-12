# Exercise 1 - reconstructing genomes from sequencing data

# load the functions you need
source("functions.R")

# we provide a function that calculates the likelihood of a certain sequence of bases; 
# this function is called calcGenoLikes and takes 4 paramaters: the sequence itself, the first allele, the secondo allele, the sequencing error rate, a boolean indicating whether the results should be returned in logarithmic scale (TRUE) or not (FALSE)

# for instance assuming that your sequence is "AATATTA", your alleles are "A" and "T", and your sequencing error rate is 0.05, then the likelihood (not in logarithms) for each genotype is given by:
# calcGenoLikes("AATATTA", "A", "T", 0.05, FALSE)

# complete all the following tasks; what it is crucial here is to not recalculate quantities that you have already computed; the idea is that you should be able to understand whether the likelihood or the prior is the same between different scenarios.

# TASK A)
# using Bayes' theorem, write the formula for the posterior probability of genotype G being AA given the sequencing data D
# write the explicit denominator assuming that your alleles are A and T and all possible genotypes are only AA, AT, TT

# TASK B)
# assuming that your data is "AAAT", your alleles are A and T, and the sequencing error rate is 0.01,
# calculate genotype posterior probability using a uniform prior:  e.g. P(G=AA) = P(G=AT) = P(G=TT) = ?

# TASK C) 
# with the same assumptions as Task B, 
# calculate genotype posterior probability using a prior based on the Hardy Weinberg Equilibrium with a frequency of T of 0.1.
# what is the Hardy Weinberg Equilibrium?
# https://en.wikipedia.org/wiki/Hardy–Weinberg_principle
# if f is the frequency of allele A, then the prior probability for all genotypes are:
# p(AA)=f^2
# p(AT)=2*f*(1-f)
# p(TT)=(1-f)^2

# do you need to calculate a new likelihood or is it the same as Task B?

# TASK D) 
# with the same assumptions as Task C, 
# calculate genotype posterior probability using a prior based on the Hardy Weinberg Equilibrium with a frequency of T of 0.1 and inbreeding coefficient of 0.2.
# what is the inbreeding coefficient?
# https://en.wikipedia.org/wiki/Inbreeding
# we need to modify our previous priors:
# if f is the frequency of allele A and I is the inbreeding coefficient, then the prior probability for all genotypes are:
# p(AA)=f^2 + I*f*(1-f)
# p(AT)=2*f*(1-f) * (1-I)
# p(TT)=(1-f)^2 + I*f*(1-f)

# do you need to calculate a new likelihood or is it the same as Tasks B-C?

# TASK E)
# with the same assumptions as Task D (same priors) but with a sequencing error rate of 0.05, calculate genotype posterior probabilities

# do you need to calculate a new likelihood or is it the same as Tasks D?

# TASK F)

# Plot all previous results (e.g. use a barplot with the 3 posterior probabilities for each scenario B-E

# TASK G) 
# What happens if we have more data? what is the confidence in our genotype inferences?
# e.g. bases <- "AAATATAAAAAAATTTTTAAATTA" and use prior as in Task B

# comment on the results...

# TASK H) 
# What happens if we have a LOT of data?
# e.g. bases <- paste(c(rep("A",1e3),rep("T",1e3)), sep="", collapse="")

# comment on the results...

# what happened? it is convenient to use numbers in log-scale: log(a*b/c) = log(a) + log(b) -log(c)
# (you can do that by selection TRUE as the fourth parameter in the calcGenoLikes function)

# what is the effect of the prior here?

# PS: if you want to calculate proper probability (in log) you have to approximate the sum of logs, see https://en.wikipedia.org/wiki/List_of_logarithmic_identities

# TASK I) optional
# how would you calculate the genotype posterior probabilities assuming that the inbreeding coefficient is not known but it can have two possible values, 0 with probability 0.7 and 0.2 with probability 0.3?
# can you write down an explicit formula for the posterior probability of genotypes AA (as in task A) under these new conditions?








