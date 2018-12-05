#######################################
# helper functions for DESeq2 analysis
#######################################

library(ggplot2)
library(gplots)
library(gtools)
library(DESeq2)
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
	rownames(results) <- res@rownames
	results$gene_id <- rownames(results)
	results$gene_name <- featureData[results$gene_id,]$gene_name
	return(results)
	}

getExpressed <- function(mat){

	     # return data frame of number of genes
	     # expressed at different sample levels

	     nsamples <- ncol(mat)
	     ngenes <- c()
	     for (i in seq(1, nsamples, 1)){
	     	 n <- nrow(mat[rowSums(mat > 0) == i,])
		 ngenes <- append(ngenes, n)
             }
	     df <- data.frame(nsamples=seq(1, nsamples, 1), ngenes=ngenes)
	     return(df)
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
	plot2 <- plot1 + geom_boxplot(outlier.alpha=0)
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

plotPCAWithCovariate <- function(pc, covariate, pcs=c("PC1", "PC2")){

        # covariate must be in same order as pc rownames

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

	# add covariate
	pca$covariate <- covariate

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

	plot1 <- ggplot(pca, aes_string(x=pc1, y=pc2, group="condition", colour="condition", shape="covariate"))
	plot2 <- plot1 + geom_point(size=3)
	plot3 <- plot2 + theme_bw() + stat_ellipse(geom="polygon", alpha=0.1, level=0.8, aes(fill=condition))
	plot4 <- plot3 + xlab(xlabel) + ylab(ylabel)
	return(plot4) 
	}


########################################
# heatmaps
########################################

heatmapMatrix <- function(mat, distfun="euclidean", clustfun="ward.D2"){

	      distf <- function(x) dist(x, method=distfun)
	      clustf <- function(x) hclust(x, method=clustfun)

	      colours <- colorRampPalette(c("blue", "white", "red"))(75)
	      mat.s <- data.frame(t(apply(mat, 1, scale)))
	      rownames(mat.s) <- rownames(mat)
	      colnames(mat.s) <- colnames(mat)
	      mat.s <- na.omit(mat.s)
	      mat.s <- mat.s[,mixedsort(colnames(mat.s))]
	      heatmap.2(as.matrix(mat.s),
	                col=colours,
			Rowv=T,
			Colv=T,
			trace="none",
			margins=c(15,15),
			hclustfun=clustf,
			distfun=distf)
	      }

heatmapMatrixWithSampleAnnotation <- function(mat, sample.annotation){

	      colours <- colorRampPalette(c("blue", "white", "red"))(75)
	      mat.s <- data.frame(t(apply(mat, 1, scale)))
	      rownames(mat.s) <- rownames(mat)
	      colnames(mat.s) <- colnames(mat)
	      mat.s <- na.omit(mat.s)
	      mat.s <- mat.s[,mixedsort(colnames(mat.s))]
	      pheatmap(as.matrix(mat.s),
	      color=colours,
	      show_rownames=F,
	      clustering_distance_cols="correlation",
	      clustering_method="complete",
	      scale="none",
	      annotation_col=sample.annotation)
	      }

################################
# MAplot
################################

getContrast <- function(dds, factor=c("condition"), contrast=c("PSC", "UC")){

	    contrast <- append(factor, contrast)
	    res <- results(dds, contrast=contrast)
	    return(res)
	    }
	    

MAPlot <- function(res, lfc=1, test.in=F, test.set, title="default title"){

         dat <- data.frame(res@listData)	    
	 rownames(dat) <- res@rownames
	 dat$gene_id <- rownames(dat)

	 if (test.in == TRUE){
	         dat$significant <- ifelse(dat$padj < 0.05 & abs(dat$log2FoldChange) > lfc & dat$gene_id %in% test.set & !(is.na(dat$padj)) & !(is.na(dat$log2FoldChange)), "Yes", "No")}
	 else{
	         dat$significant <- ifelse(dat$padj < 0.05 & abs(dat$log2FoldChange) > lfc & !(is.na(dat$padj)) & !(is.na(dat$log2FoldChange)), "Yes", "No")}

       nup <- nrow(dat[dat$significant=="Yes" & dat$log2FoldChange > lfc & !(is.na(dat$padj)) & !(is.na(dat$log2FoldChange)),])
       ndown <- nrow(dat[dat$significant=="Yes" & dat$log2FoldChange < (-lfc) & !(is.na(dat$padj)) & !(is.na(dat$log2FoldChange)),])

       nup <- paste("Upregulated = ", nup, sep="")
       ndown <- paste("Downregulated = ", ndown, sep="")

       plot1 <- ggplot(na.omit(dat), aes(x=log2(baseMean + 1), y=log2FoldChange, colour=significant))
       plot2 <- plot1 + geom_point(pch=18, size=1)
       plot3 <- plot2 + scale_colour_manual(values=c("grey", "darkRed"))
       plot4 <- plot3 + theme_bw()
       plot5 <- plot4 + theme(text=element_text(size=10))
       plot6 <- plot5 + geom_hline(yintercept=c(-1,1,0), linetype="dashed")
       plot7 <- plot6 + xlab("Mean expression") + ylab("Log2 fold change")
       plot8 <- plot7 + annotate("text", x=9, y=6, label=nup, size=3)
       plot9 <- plot8 + annotate("text", x=9, y=-6, label=ndown, size=3)
       plot10 <- plot9 + ggtitle(title)
       return (plot10)
       }


#########################################
#########################################
# transformations
#########################################
#########################################

log2cpm <- function(counts){

	cpm <- log2(sweep(counts, 2, colSums(counts)/1000000, "/") + 1)
	return(cpm)
	}

relab <- function(counts){

      relab <- (sweep(counts, 2, colSums(counts), "/"))*100
      return(relab)
      }