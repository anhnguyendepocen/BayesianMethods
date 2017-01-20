# Solutions exercise 1 - reconstructing genomes from sequencing data

source("functions.R")

bases <- "AAAT"

# calculate genotype posterior probability with different priors:

# P(G=AA|D) = P(D|G=AA) * P(G=AA) / P(D)

# A) uniform prior:  P(G=AA) = P(G=AT) = P(G=TT) = 1/3

likes <- calcGenoLikes("AAAT", "A", "T", 0.01, FALSE)
priors_A <- rep(1/3, 3)
posteriors_A <- likes * priors_A / sum(likes * priors_A)

# B) HWE-based prior assuming a frequency of T of 0.3
f <- 0.1
priors_B <- c( (1-f)^2, 2*f*(1-f), f^2 )
posteriors_B <- likes * priors_B / sum(likes * priors_B)

# C) Inbreeding of I=0.2
I <- 0.2
priors_C <- c( (1-f)^2 + I*f*(1-f), 2*f*(1-f) * (1-I), f^2 + I*f*(1-f))
posteriors_C <- likes * priors_C / sum(likes * priors_C)

# D) Inbreeding and higher error rate 0.05
likes_D <- calcGenoLikes("AAAT", "A", "T", 0.05, FALSE)
posteriors_D <- likes_D * priors_C / sum(likes_D * priors_C)

# plot
posteriors <- matrix(c(posteriors_A, posteriors_B, posteriors_C, posteriors_D), ncol=3, nrow=4, byrow=T)
barplot(posteriors, beside=T, legend.text=toupper(letters[1:4]), names=c("AA","AT","TT"), ylab="Posterior probability", xlab="Genotype")

# E) what if we have more data? what is the confidence here?

bases <- "AAATATAAAAAAATTTTTAAATTA"
likes_E <- calcGenoLikes(bases, "A", "T", 0.05, FALSE)
posteriors_E <- likes_E * priors_C / sum(likes_E * priors_C)

# F) what if we have a LOT of data?

bases <- paste(c(rep("A",1e3),rep("T",1e3)), sep="", collapse="")
likes_F <- calcGenoLikes(bases, "A", "T", 0.05, FALSE)
posteriors_F <- likes_F * priors_C / sum(likes_F * priors_C)

# what happened? it is convenient to use numbers in log-scale: log(a*b/c) = log(a) + log(b) -log(c)

loglikes_F <- calcGenoLikes(bases, "A", "T", 0.05, TRUE)
loglikes_F + log(priors_C)

# what is the effect of the prior here?

# PS: if you want to calculate proper probability (in log) you have to approximate the sum of logs, see https://en.wikipedia.org/wiki/List_of_logarithmic_identities


