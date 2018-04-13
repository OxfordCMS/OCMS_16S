################################################
################################################
################################################
# helper functions for WGCNA analysis
################################################
################################################
################################################


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

	 
sampleMatrix <- function(cor.matrix, n, nsamples=1000){

	     mean.cors = c()
	     # sample a correlation matrix
	     for (i in 1:nsamples){
	     	 ids <- sample(rownames(cor.matrix), n, replace=FALSE)
		 cor.sub <- cor.matrix[ids, ids]
		 diag(cor.sub) <- NA
		 mean.cor <- mean(cor.sub, na.rm=T)
		 mean.cors <- append(mean.cors, mean.cor)
	     }
	     return(mean.cors)
	     }
