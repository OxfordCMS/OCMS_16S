##############################################
# functions for dealing with things related to
# GO.py
##############################################

library(data.table)

getGenesetGenes <- function(genesets, geneset){

		genesets <- genesets[genesets$V3 == geneset,]
		genes <- genesets$V2
		return(genes)
		}

plotPathways <- function(pathway.results){

	     p.dat <- pathway.results
	     p.dat <- p.dat[p.dat$code=="+",]
	     p.dat <- p.dat[order(p.dat$ratio, decreasing=T),]
	     p.dat$label <- paste(round(p.dat$fdr, 3), " (", sep="")
	     p.dat$label <- paste(p.dat$label, p.dat$scount, sep="")
	     p.dat$label <- paste(p.dat$label, ")", sep="")

	     p.dat$goid <- factor(p.dat$goid, levels=p.dat$goid)

	     plot10 <- ggplot(p.dat, aes(x=goid, y=ratio, label=scount))
	     plot11 <- plot10 + geom_bar(stat="identity") + geom_text(vjust=-1)
	     plot12 <- plot11 + theme_bw() + theme(text=element_text(angle=90))
	     plot13 <- plot12 + ylim(c(0,5))
	     plot14 <- plot13 + xlab("") + ylab("Fold enrichment")
	     return(plot14)
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