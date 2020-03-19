####################################################
# basic summary plots of data coming from
# readqc pipeline - access via sqlite database
####################################################

library("RSQLite")
library("ggplot2")
library("reshape")

##########
# readQC
##########

getSamples <- function(db){
	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      samples <- dbGetQuery(db, 'SELECT track FROM status_summary as s WHERE track NOT LIKE "%%trimmed%%"')
              rsamples = c()
              for (s in samples){
                  rsamples <- append(gsub("-", "_",s), rsamples)
		  }
	      return(rsamples)
	      }
getSampleQC <- function(db, sample, metric){
	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      toget = paste(sample, "fastqc", sep="_")
	      toget = paste(toget, metric, sep="_")
	      df <- dbGetQuery(db, paste('SELECT * FROM ', toget, sep=""))
	      df$condition <- unlist(strsplit(sample, "_"))[2]
	      df$track <- sample
	      return(df)
	      }

combineSampleQC <- function(db, metric, samples){
		datalist = list()
		for (i in 1:length(samples)){
		    df <- getSampleQC(db, samples[i], metric)
		    datalist[[i]] <- df
		}
		combined <- do.call(rbind, datalist)
		return(combined)
		}

plotPerBaseSequenceQuality <- function(df){
			   plot1 <- ggplot(df, aes(x=Base,
			                           y=Mean,
						   group=track,
						   colour=condition))
			   plot2 <- plot1 + geom_line() 
			   plot3 <- plot2 + theme_bw()
			   plot4 <- plot3 + ylim(c(0, 40)) + xlim(1, max(df$Base))
			   plot5 <- plot4 + ggtitle("Per base squence quality")
			   return(plot5)
			   }

plotPerSequenceQuality <- function(df){
                      plot1 <- ggplot(df, aes(x=Quality,
		                              y=Count,
					      group=track,
					      colour=condition))
		      plot2 <- plot1 + geom_line()
		      plot3 <- plot2 + theme_bw()
		      plot4 <- plot3 + ggtitle("Per base squence quality")
		      return(plot4)
		      }

plotPerSequenceGC <- function(df){
                      plot1 <- ggplot(df, aes(x=GC_Content,
		                              y=Count,
					      group=track,
					      colour=condition))
		      plot2 <- plot1 + geom_line()
		      plot3 <- plot2 + theme_bw()
		      plot4 <- plot3 + ggtitle("Per sequence GC content")
		      return(plot4)
		      }

plotSequenceDuplicationLevels <- function(df){
			      plot1 <- ggplot(df, aes(x=Duplication_Level,
			                              y=Percentage_of_total,
			                              group=track,
						      colour=condition))
			      plot2 <- plot1 + geom_line()
			      plot3 <- plot2 + theme_bw()
			      plot4 <- plot3 + ggtitle("Per sequence duplication levels")
			      return(plot4)
			      }

##########
# Mapping
##########

getReadCounts <- function(db){
	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      df <- dbGetQuery(db, 'SELECT * FROM reads_summary;')
	      return(df)
	      }
	      
plotReadsSummary <- function(df){

	plot1 <- ggplot(df, aes(x=track, y=total_reads))
	plot2 <- plot1 + geom_bar(stat="identity") + theme_bw()
	plot3 <- plot2 + theme(axis.text.x=element_text(angle=90))
	plot4 <- plot3 + ggtitle("Total reads")
	return(plot4) 
	}

getConds <- function(samples){

	 conds <- unlist(strsplit(samples, "-"))
	 conds <- conds[seq(2, length(conds), 3)]
	 return(conds)
	 }

plotReadsSummaryByCondition <- function(df){

			    df$cond <- getConds(df$track)
			    plot1 <- ggplot(df, aes(x=cond, y=total_reads, group=cond))
			    plot2 <- plot1 + geom_boxplot() + theme_bw()
			    plot3 <- plot2 + theme(axis.text.x=element_text(angle=90))
			    plot4 <- plot3 + ggtitle("Read count by condition")
			    return(plot4)
			    }

getProportionReadsMapped <- function(db){

	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      df <- dbGetQuery(db, 'SELECT track, reads_total, reads_unmapped FROM bam_stats;')
	      df$p <- ((df$reads_total - df$reads_unmapped)/df$reads_total)*100
	      return(df)
	      }

plotProportionReadsMapped <- function(df){

			  plot1 <- ggplot(df, aes(x=track, y=p))
			  plot2 <- plot1 + geom_bar(stat="identity")
			  plot3 <- plot2 + theme_bw()
			  plot4 <- plot3 + theme(axis.text.x=element_text(angle=90))
			  plot5 <- plot4 + ggtitle("Proportion of reads mapped")
			  return(plot5)
			  }

getNumberSecondaryAlignments <- function(db){

	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      df <- dbGetQuery(db, 'SELECT track, alignments_mapped, alignments_secondary FROM bam_stats;')
	      df <- melt(df)
	      dbDisconnect(db)
	      return(df)
	      }

plotNumberSecondaryAlignments <- function(df){

			      plot1 <- ggplot(df, aes(x=track, y=value, fill=variable))
			      plot2 <- plot1 + geom_bar(stat="identity")
			      plot3 <- plot2 + theme_bw()
			      plot4 <- plot3 + scale_fill_manual(values=c("grey", "red3"))
			      plot5 <- plot4 + ggtitle("Secondary alignments")
			      plot6 <- plot5 + theme(axis.text.x=element_text(angle=90))
			      return(plot6)
			      }
			      
getContextStats <- function(db){

	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      df <- dbGetQuery(db, 'SELECT track, repeats, intergenic, intron, three_prime_utr, five_prime_utr, protein_coding, total FROM context_stats;')

	      # make counts into proportions (by row)
	      df2 <- df[ ,grep("track", colnames(df), invert=T)]

	      # remove total count
	      df2 <- df2[ ,grep("total", colnames(df2), invert=T)]
	      df2 <- sweep(df2, 1, rowSums(df2), "/")*100
	      df2$track <- df$track
	      return(df2)
	      }

getContextCoverage <- function(db, lengths){

	      sqlite <- dbDriver("SQLite")
	      db <- dbConnect(sqlite, db)
	      df <- dbGetQuery(db, 'SELECT track, repeats, intergenic, intron, three_prime_utr, five_prime_utr, protein_coding, total FROM context_stats;')

	      df2 <- df[ ,grep("track", colnames(df), invert=T)]

	      # remove total count
	      df2 <- df2[ ,grep("total", colnames(df2), invert=T)]

	      # read lengths
	      lengths <- read.csv(lengths, header=T, stringsAsFactors=F, sep="\t", row.names=1)

	      # make columns of df2 the same as rows of lengths
	      df2 <- df2[,rownames(lengths)]

	      # reads per million
	      rpm <- data.frame(sweep(df2, 1, (rowSums(df2)/1000000), "/"))

	      # reads per kilobase
	      rpkm <- data.frame(t(apply(rpm, 1, function(x) x/(lengths$length/1000))))

	      # add track back
	      rpkm$track <- df$track

	      return(rpkm)
	      }	      
		   

plotContextStats <- function(df){

		rownames(df) <- df$track
		df <- df[,grep("track", colnames(df), invert=T)]	

		df$track <- rownames(df)

		# set colours
		cols <- rainbow(ncol(df)-1, v=0.8, s=0.5)

		# reshape data frame
		df.m <- melt(df)

		plot1 <- ggplot(df.m, aes(x=track, y=value))
		plot2 <- plot1 + geom_bar(stat="identity")
		plot3 <- plot2 + theme_bw()
		plot4 <- plot3 + scale_fill_manual(values=cols)
		plot5 <- plot4 + theme(axis.text.x=element_text(angle=90))
		return(plot5 + facet_wrap(~variable))
		}
		 