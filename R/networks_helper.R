################################################
################################################
################################################
# helper functions for co-expression
# network analysis
################################################
################################################
################################################

readList <- function(infile){

	 dat <- read.csv(infile, header=F, sep="\t", stringsAsFactors=F)
	 g <- dat$V1
	 return(g)
	 }

log2cpm <- function(counts){

	cpm <- log2(sweep(counts, 2, colSums(counts)/1000000, "/") + 1)
	return(cpm)
	}

residualise <- function(row, conds, scale=T){

	    # save standardised residuals after regressing
	    # out effect of conds

	    linmod <- lm(unlist(row)~as.factor(conds))
	    residuals <- residuals(linmod)
	    if (scale==TRUE){
	       residuals <- scale(residuals)}
	    else{
		residuals <- residuals
	    }
	    return(residuals)
	    }

	 
sampleMatrix <- function(cor.matrix, n, fun="mean", nsamples=1000){

	     mean.cors = c()
	     # sample a correlation matrix
	     for (i in 1:nsamples){
	     	 ids <- sample(rownames(cor.matrix), n, replace=FALSE)
		 cor.sub <- cor.matrix[ids, ids]
		 diag(cor.sub) <- NA
		 if (fun == "mean"){
		     mean.cor <- mean(cor.sub, na.rm=T)}
		 else{
		     mean.cor <- median(cor.sub, na.rm=T)
		 }
             mean.cors <- append(mean.cors, mean.cor)
	     }
	     return(mean.cors)
	     }


calcP <- function(expected, observed){

      # expected = is a vector of means from
      # sampling a network
      # observed = single mean value from subnetwork
      # of interest
      
      above <- ifelse(expected >= observed, TRUE, FALSE)
      nabove <- length(above[above == TRUE])
      nexpected <- length(expected)
      p <- nabove/nexpected
      return(p)
      }
