#######################################
# helper functions for DESeq2 analysis
#######################################

library(ggplot2)
library(gplots)
library(gtools)
library(RSQLite)

#################################
# making conditions and features
#################################

makeConds <- function(samples){
	  
	    conditions <- unlist(strsplit(samples, "\\."))
	    conditions <- conditions[seq(2,length(conditions),3)]
	    return(conditions)
	    }
	    
makeColData <- function(mat){

	    samples <- colnames(mat)
	    conditions <- makeConds(samples)
	    coldata <- data.frame(sample=samples, condition=conditions)
	    return(coldata)
	    }

makeFeatureData <- function(db, mat){

		sqlite <- dbDriver("SQLite")
	      	db <- dbConnect(sqlite, db)
	      	ens2gene <- dbGetQuery(db, 'SELECT distinct(gene_id), gene_name FROM gene_info')
		rownames(ens2gene) <- ens2gene$gene_id
		ens <- rownames(mat)
		genes <- ens2gene[ens,]$gene_name
		dbDisconnect(db)
		featuredata <- data.frame(gene_id = ens, gene_name = genes)
		return (featuredata)
		}

getResultsTable <- function(res, featureData){

	rownames(featureData) <- featureData$gene_id
	results <- data.frame(res@listData)
	results$gene_id <- rownames(results)
	results$gene_name <- featureData[results$gene_id,]$gene_name
	return(results)
	}

plotMeanSd <- function(mat){

		 means <- rowMeans(mat)
		 sds <- apply(mat, 1, sd)
		 plot(means, sds, col="grey", pch=16)

		 # smooth line
		 s <- smooth.spline(means, sds)
		 lines(s, col="red", lty=2)
		 }

plotGOI <- function(mat, goi=""){

	# some funky behaviour with "-"
	rownames(mat) <- gsub("-", "_", rownames(mat))
	goi <- gsub("-", "_", goi)

	conds <- makeConds(colnames(mat))
	dat <- mat[goi,]
	df <- data.frame(t(dat))
	df$cond <- conds
	plot1 <- ggplot(df, aes_string(x="cond", y=goi, colour="cond"))
	plot2 <- plot1 + geom_boxplot()
	plot3 <- plot2 + geom_jitter()
	plot4 <- plot3 + theme_bw()
	plot5 <- plot4 + ggtitle(goi)
	return(plot5)
}
	

#################################
# principle components analysis
#################################

PCA <- function(df){

    pc <- prcomp(t(df))
    return (pc)
    }

getPCA <- function(pc){

       pcs <- data.frame(pc$x)
       return(pcs)
       }

getVE <- function(pc, component="PC1"){

      pve <- summary(pc)$importance[,component][[2]]
      return (pve)
      }

plotCumulativeProportion <- function(pc){

			 cp <- data.frame(summary(pc)$importance["Cumulative Proportion",])
			 colnames(cp) <- "cp"
			 cp$PC <- rownames(cp)
			 cp$PC <- factor(cp$PC, levels=cp$PC)
			 plot1 <- ggplot(cp, aes(x=PC, y=cp))
			 plot2 <- plot1 + geom_bar(stat="identity")
			 plot3 <- plot2 + theme_bw()
			 plot4 <- plot3 + ggtitle("Cumulative Proportion")
			 return (plot4)
			 }

plotPCA <- function(pc, pcs=c("PC1", "PC2")){

	# get variance explained for each component
	ve1 <- getVE(pc, component=pcs[1])
	ve2 <- getVE(pc, component=pcs[2])

	ve1 <- round(ve1, 2)*100
	ve2 <- round(ve2, 2)*100

	# get data frame of components
	pca <- data.frame(pc$x)

	# add conditions
	pca$condition <- makeConds(rownames(pca))
	pca$condition <- factor(pca$condition, levels=unique(pca$condition))

	# plot
	
	pc1 <- pcs[1]
	pc2 <- pcs[2]

	# labels
	xlabel <- paste(pc1, ve1, sep=" (")
	xlabel <- paste(xlabel, "%", sep="")
	xlabel <- paste(xlabel, ")", sep="")
	ylabel <- paste(pc2, ve2, sep=" (")
	ylabel <- paste(ylabel, "%", sep="")	
	ylabel <- paste(ylabel, ")", sep="")

	plot1 <- ggplot(pca, aes_string(x=pc1, y=pc2, group="condition", colour="condition"))
	plot2 <- plot1 + geom_point(size=3)
	plot3 <- plot2 + theme_bw() + stat_ellipse(geom="polygon", alpha=0.1, level=0.8, aes(fill=condition))
	plot4 <- plot3 + xlab(xlabel) + ylab(ylabel)
	return(plot4) 
	}

########################################
# heatmap
########################################

heatmapMatrix <- function(mat){

	      colours <- colorRampPalette(c("blue", "white", "red"))(75)
	      mat.s <- data.frame(t(apply(mat, 1, scale)))
	      rownames(mat.s) <- rownames(mat)
	      colnames(mat.s) <- colnames(mat)
	      mat.s <- na.omit(mat.s)
	      mat.s <- mat.s[,mixedsort(colnames(mat.s))]
	      heatmap.2(as.matrix(mat.s), col=colours, Rowv=T, Colv=T, trace="none", margins=c(15,15))
	      }

