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

