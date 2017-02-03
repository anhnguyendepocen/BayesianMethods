## R functions

library("fields")


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


# plot 2D SFS
plot2DSFS<-function(sfs, xlab="", ylab="", main="") {

	# this function plots a 2D-SFS in log10 scale; it requires "fields" package
	nchroms_pop1<-nrow(sfs)-1
        nchroms_pop2<-ncol(sfs)-1

	sfs=log10(sfs)
	sfs[which(is.na(sfs))]=0
	sfs[which(sfs==(-Inf))]=0	

	brk=seq(min(sfs),max(sfs)+0.3,0.3)
	lab.brk=round(10^(seq(min(sfs),max(sfs)+0.3,0.3)))
	lab.brk[1]=0
	cols=rainbow(length(brk)-1)
	cols[1]="white"

	image.plot(x=seq(0,nchroms_pop1), y=seq(0,nchroms_pop2), z=sfs, breaks=brk, col=cols, lab.breaks=lab.brk, xlab=xlab, ylab=ylab, main=main)


}


# fold 2D SFS
fold2DSFS<-function(unfolded) {

	# this function takes as input an unfolded spcetrum and returns a folded one

        folded<-unfolded

        sample_size_pop1<-(nrow(unfolded)-1)/2
        sample_size_pop2<-(ncol(unfolded)-1)/2

        freq_at50<-sample_size_pop1+sample_size_pop2
 
        for (i in 1:nrow(unfolded)) {
                for (j in 1:ncol(unfolded)) {
                        daf1=(i-1) ##take the values chosen from the rows and subtract 1 for the true daf
                        daf2=(j-1) ##take the values chosen from the columns and subtract 1 for the true daf
                        if ((daf1+daf2)>freq_at50) { ## when value is major
                                minor1=(2*sample_size_pop1)-daf1; ## converts to minor
                                minor2=(2*sample_size_pop2)-daf2; ## converts to minor
                                folded[(minor1+1),(minor2+1)]=folded[(minor1+1), (minor2+1)]+unfolded[(daf1+1), (daf2+1)]; ##add 1 back to the values to get to the true spot in the matrix and replace it by adding the values that were in the majors to the values in the minors
                                folded[(daf1+1),(daf2+1)]=NA ##replace all values where majors originally were with NA
                        }
                }
        }

        folded

}


# calculate FST
doFST<-function(sfs) {

	# this function compute the FST from a 2D-SFS, using the Reynold's et al estimator

        sfs=sfs/sum(sfs, na.rm=T)
        nind1=(nrow(sfs)-1)/2
        nind2=(ncol(sfs)-1)/2

        nums=denoms=fsts=matrix(NA, ncol=ncol(sfs), nrow=nrow(sfs))
        for (i in 1:nrow(sfs)) {
                for (j in 1:ncol(sfs)) {
                        f1=(i-1)/(nind1*2); f2=(j-1)/(nind2*2)
                        tmp=reynolds(f1,f2,nind1*2,nind2*2)
                        nums[i,j]=tmp[1]
                        denoms[i,j]=tmp[2]
                        fsts[i,j]=tmp[3]
                }
        }

        fst=sum(sfs*nums,na.rm=T)/sum(sfs*denoms,na.rm=T)
        fst
}


# FST estimator
reynolds<-function(pl1,pl2,n1,n2) {

	# this functions compute the FST estimator from Reynolds et al. from sample allele frequencies
	somma=sommaden=0;

	alfa1=1-((pl1^2)+((1-pl1)^2))
	alfa2=1-((pl2^2)+((1-pl2)^2))
	Al = (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) - (((n1+n2)*(n1*alfa1+n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))
	AlBl= (0.5*(((pl1-pl2)^2)+(((1-pl1)-(1-pl2))^2))) + (((4*n1*n2 - n1 - n2)*(n1*alfa1 + n2*alfa2)) / ((4*n1*n2)*(n1+n2-1)))

	if (!is.na(Al) & !is.na(AlBl)) {
  		somma=somma+Al
  		sommaden=sommaden+AlBl
 	}

	if (somma==0 & sommaden==0) {
      		reyn=NA
 	} else {
        	reyn=somma/sommaden
        	if(reyn<0) reyn=0
 	}

 	c(somma, sommaden, reyn)

}


fromMStoSFSwith1site<-function(msfile, nr_repetitions, nr_chroms_pop1, nr_chroms_pop2 ) {

	# this functions read ms output file (if only 1 site is simulated) and returns a matrix containing the 2D-SFS
	# originally written by S.D. Gopal

        nr_chromosomes<-nr_chroms_pop1+nr_chroms_pop2

	# read ms file
        sequencedata<-readLines(msfile)##read the ms command output into R
        sequencedata<-suppressWarnings(as.numeric(sequencedata)) ##make output numeric
        sequencedata<-sequencedata[!is.na(sequencedata)] ##remove NAs which are words in this case

        sequencedata_matrix<-matrix(sequencedata, nrow= (nr_chromosomes), ncol= (nr_repetitions))

        chroms_pop1<-sequencedata_matrix[1:nr_chroms_pop1, (1:nr_repetitions)] ### subset pop1 data
        chroms_pop2<-sequencedata_matrix[(nr_chroms_pop1+1):(nr_chroms_pop1+nr_chroms_pop2), (1:nr_repetitions)] ### subset pop2 data

        sfs<-matrix(0, nrow= (nr_chroms_pop1+1), ncol= (nr_chroms_pop2+1)) ### create a matrix filled with zeroes

        for (j in 1:nr_repetitions) {

                ### sum the value across each site in pop1 and pop2 respectively to calculate the derived allele frequencies 

                daf_pop1<-sum(chroms_pop1[,j]) + 1 # +1 to convert into coordinates
                daf_pop2<-sum(chroms_pop2[,j]) + 1

                sfs[daf_pop1, daf_pop2]<-sfs[daf_pop1, daf_pop2]+1 ## place summations in the matrix 

        }

        sfs


}






