## R functions

# calculate genotype likelihoods
calcGenoLikes <- function(bases, major, minor, errorRate=0.01, log.scale=TRUE) {

        # initialise
        if (log.scale) likes <- rep(0, 3) else likes <- rep(1, 3)

        alleles <- c('A', 'C', 'G', 'T')
        iter <- 0
        reads <- unlist(strsplit(bases, split=""))

        # cycle across possible genotypes
        for ( j in list( c(match(major, alleles),match(major, alleles)), c(match(major, alleles),match(minor, alleles)), c(match(minor, alleles),match(minor, alleles)) )  ) {

                iter <- iter + 1

                # cycle across all reads
                for (i in  1:length(reads)) {

                        sublike = 0.0

                        if (alleles[j[1]] == reads[i]) sublike = sublike + (1-errorRate)/2 else sublike = sublike + (errorRate/3)/2

                        if (alleles[j[2]] == reads[i]) sublike = sublike + (1-errorRate)/2 else sublike = sublike + (errorRate/3)/2

                        if (log.scale) likes[iter] = likes[iter] + log(sublike)  else  likes[iter] = likes[iter] * sublike

                }
        }

        names(likes) <- c( paste(c(major,major), sep="", collapse=""), paste(c(major,minor), sep="", collapse=""), paste(c(minor,minor), sep="", collapse="") )
        likes
}


