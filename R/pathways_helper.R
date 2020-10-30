##############################################
# functions for dealing with things related to
# GO.py
##############################################

library(data.table)
library(pheatmap)

getGenesetGenes <- function(genesets, geneset){

		genesets <- genesets[genesets$V3 == geneset,]
		genes <- genesets$V2
		return(genes)
		}

plotPathways <- function(pathway.results, plot.all=FALSE, relevel=TRUE, text.size=6, x.label="description"){

	     p.dat <- pathway.results
	     if (plot.all == TRUE){
	       p.dat <- p.dat
	     }
	     else
	     {
	     p.dat <- p.dat[p.dat$code=="+",]
	     }
	     p.dat <- p.dat[order(p.dat$ratio, decreasing=T),]
	     p.dat$label <- paste(round(p.dat$fdr, 3), " (", sep="")
	     p.dat$label <- paste(p.dat$label, p.dat$scount, sep="")
	     p.dat$label <- paste(p.dat$label, ")", sep="")

	     if (relevel == TRUE){
	        p.dat$goid <- factor(p.dat$goid, levels=p.dat$goid)
	     }
	     plot1 <- ggplot(p.dat, aes(x=goid, y=ratio, label=scount, fill=goid)) 
	     plot2 <- plot1 + geom_bar(stat="identity") + geom_text(vjust=-1, size=text.size)
	     plot3 <- plot2 + theme_bw() + theme(axis.text.x=element_text(angle=90))
	     plot4 <- plot3 + xlab("") + ylab("Fold enrichment") + scale_fill_manual(values=rep("grey", nrow(p.dat))) + theme(legend.position = "none")
	     return(plot4)
	     }

buildGenesetGenelist <- function(pathway.table, genesets){

		genes <- c()
		for (i in 1:nrow(pathway.table)){
		    geneset <- pathway.table$goid[i]
		    gs <- as.character(getGenesetGenes(genesets, geneset))
		    genes <- append(gs, genes)
		}
		genes <- unique(genes)
		return(genes)
		}

buildGeneMatrix <- function(total, pathway.table, genesets){

		df.list <- list()
		for (i in 1:nrow(pathway.table)){
		    geneset <- pathway.table$goid[i]
		    genes <- getGenesetGenes(genesets, geneset)
		    gs <- c()
		    values <- c()
		    for (g in total){
		        gs <- append(g, gs)
			if (g %in% genes){
			   values <- append(1, values)}
			else{
			   values <- append(0, values)}
	            }
		    v <-  geneset
		    df <- data.frame(gene=gs, value=values)
		    colnames(df) <- c("gene", v)
		    df.list[[i]] <- data.table(df)
		    
		 }
		 full.mat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by="gene", all=T), df.list)
		 full.mat <- as.data.frame(full.mat)
		 rownames(full.mat) <- full.mat$gene
		 full.mat <- full.mat[,2:ncol(full.mat)]
		 return(full.mat)
}


binaryHeatmap <- function(mat, colour="blue", labels=TRUE){

    # draw heatmap

    cols <- colorRampPalette(c("white", colour))(2)
    distfun=function(x) dist(x, method="binary")
    pheatmap(mat, color=cols, cluster_distance_rows=distfun, cluster_distance_cols=distfun, show_rownames = labels)
}

