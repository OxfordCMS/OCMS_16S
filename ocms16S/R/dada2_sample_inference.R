###################################################
###################################################
# Error model, de-replication and sample inference
###################################################
###################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dada2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("futile.logger"))

# make options list
option_list <- list(
               make_option(c("--filtF"), default=NA, type="character",
                           help="filtered forward fastq file [default %default]"),
               make_option(c("--filtR"), default=NA, type="character",
                           help="filtered reverse fastq file [default %default]"),
	       make_option(c("-n", "--nbases"), default=10000000,
                           help="number of reads to learn error model [default %default]"),
	       make_option(c("--omega-a"), default=1e-40,
                           help="p-value for new clusters [default %default]"),
	       make_option(c("--use-quals"), default=TRUE,
                           help="If TRUE, the dada(...) error model takes into account
                           the consensus quality score of the dereplicated unique sequences.
                           If FALSE, quality scores are ignored. Default is TRUE.
                           [default %default]"),
               make_option(c("--use-kmers"), default=TRUE,
                           help="If TRUE, a 5-mer distance screen is performed prior to
                           performing each pairwise alignment, and if the 5mer-distance is
                           greater than KDIST_CUTOFF, no alignment is performed. Default is
                           TRUE. [default %default]"),
	       make_option(c("--kdist-cutoff"), default=0.42,
                           help="The default value of 0.42 was chosen to screen pairs
                           of sequences that differ by >10%, and was calibrated on Illumina
                           sequenced 16S amplicon data. The assumption is that sequences that
                           differ by such a large amount cannot be linked by amplicon errors
                           (i.e. if you sequence one, you won't get a read of other) and so
                           careful (and costly) alignment is unnecessary. [default %default]"),
	       make_option(c("--band-size"), default=16,
                           help="When set, banded Needleman-Wunsch alignments are
                           performed. Banding restricts the net cumulative number of
                           insertion of one sequence relative to the other. The default value
                           of BAND_SIZE is 16. If DADA is applied to marker genes with high
                           rates of indels, such as the ITS region in fungi, the BAND_SIZE
                           parameter should be increased. Setting BAND_SIZE to a negative
                           number turns off banding (i.e. full Needleman-Wunsch). [default %default]"),
	       make_option(c("--gap-penalty"), default=-8,
                           help="The cost of gaps in the Needleman-Wunsch alignment.
                           Default is -8.[default %default]"),
	       make_option(c("--homopolymer-gap-penalty"), default=NULL,
                           help="The cost of gaps in homopolymer regions
                          (>=3 repeated bases). Default is NULL, which causes homopolymer
                          gaps to be treated as normal gaps. [default %default]"),
	       make_option(c("--min-fold"), default=1,
                           help="The minimum fold-overabundance for sequences to form new
                           clusters. Default value is 1, which means this criteria is
                           ignored.[default %default]"),
               make_option(c("--min-hamming"), default=1,
                           help="minimum hamming distance to define new cluster [default %default]"),
               make_option(c("--min-abundance"), default=1,
                           help="The minimum abundance for unique sequences form new
                           clusters. Default value is 1, which means this criteria is
                           ignored.[default %default]"),
               make_option(c("--max-clust"), default=0,
                           help="The maximum number of clusters. Once this many clusters
                           have been created, the algorithm terminates regardless of whether
                           the statistical model suggests more real sequence variants exist.
                           If set to 0 this argument is ignored. Default value is 0.[default %default]"),
               make_option(c("--max-consist"), default=10,
                           help="The maximum number of steps when selfConsist=TRUE. If
                           convergence is not reached in MAX_CONSIST steps, the algorithm
                           will terminate with a warning message. Default value is 10.
                           [default %default]"),
               make_option(c("--filter"), default=0,
                           help="filter sequences out if (n0+n1)/abundance < --filter. This hopes to remove
                           spurious sequences that are a partition that contains sequences with few true
                           matches [default %default]"),
               make_option(c("-o", "--outdir"), default=".",
                           help="output directory - name will be based on file name [default %default]")
			   )

# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################

opt <- parse_args(OptionParser(option_list=option_list))

# get sample name
sample.name <- gsub(".fastq.1.gz", "", basename(opt$`filtF`))

# check that filtF is specified
if (is.na(opt$`filtF`)){
   stop("--filtF parameter must be provided. See script usage (--help)")
   }

directory <- dirname(opt$`filtF`)

###########################################################################
# learn errors and do sample inference - separates out for paired-end or
# single end at the moment although this could probably be refactored in
# the future
###########################################################################

flog.info(gsub("nbases", as.character(opt$`nbases`), "learning errors on the first nbases bases"))

if (is.na(opt$`filtR`)){

   # error model inference
   errF <- learnErrors(opt$`filtF`, multithread=FALSE, nbases=opt$`nbases`)
   flog.info("plotting error model")
   p <- plotErrors(errF, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errF.png", sep="_")
   ggsave(filename, type="cairo")

   # dereplication
   flog.info("dereplication")
   derepF <- derepFastq(opt$`filtF`, verbose=TRUE)

   # sample inference
   flog.info("sample inference")	
   dadaF <- dada(derepF,
                 err=errF,
		 multithread=FALSE,
		 OMEGA_A=opt$`omega-a`,
		 USE_QUALS=opt$`use-quals`,
		 USE_KMERS=opt$`use-kmers`,
		 KDIST_CUTOFF=opt$`kdist-cutoff`,
		 BAND_SIZE=opt$`band-size`,
		 GAP_PENALTY=opt$`gap-penalty`,
		 HOMOPOLYMER_GAP_PENALTY=opt$`homopolymer-ga-penalty`,
		 MIN_FOLD=opt$`min-fold`,
                 MIN_HAMMING=opt$`min-hamming`,
                 MIN_ABUNDANCE=opt$`min-abundance`,
		 MAX_CLUST=opt$`max-clust`,
		 MAX_CONSIST=opt$`max-consist`
                 )

   # return a dataframe - I find this more intuitive
   dadaF.df <- as.data.frame(dadaF$denoised)
   dadaF.df$sequence <- rownames(dadaF.df)
   dadaF.df <- data.frame(sequence=dadaF.df$sequence, abundance=dadaF.df[,1])

   # get diagnostic clustering data frame
   df.clustering <- dadaF$clustering
   clustering.filename <- paste0(opt$`outdir`, "/", sample.name, "_clustering.tsv")
   write.table(df.clustering, clustering.filename, sep="\t", row.names=FALSE)

   # filtering options
   sequences.to.keep <- df.clustering$sequence[(df.clustering$n0 + df.clustering$n1)/df.clustering$abundance >= opt$`filter`]
   dadaF.df <- dadaF.df[dadaF.df$sequence %in% sequences.to.keep,]

   # remove chimeras
   flog.info("removing chimeric sequences")
   dadaF.df.nochim <- removeBimeraDenovo(dadaF.df, method="consensus")

   # write out table
   flog.info("writing results")
   outfile <- paste(opt$`outdir`, sample.name, sep="/")
   outfile <- paste(outfile, "_seq_abundance.tsv", sep="")
   write.table(dadaF.df.nochim, file=outfile, sep="\t", quote=F, row.names=F)

   # write summaries
   # get summaries of reads passing each stage
   flog.info("writing summaries")
   getN <- function(x) sum(getUniques(x))
   track <- data.frame(cbind(getN(dadaF), getN(dadaF.df.nochim)))
   colnames(track) <- c("denoisedF", "nonchim")
   track$sample <- sample.name

   summary.outfile <- paste(opt$`outdir`, sample.name, sep="/")
   summary.outfile <- paste(summary.outfile, "_summary.tsv", sep="")
   write.table(track, summary.outfile, sep="\t", row.names=F)


} else {

   # error model learning
   errR <- learnErrors(opt$`filtR`, multithread=FALSE, nbases=opt$`nbases`)
   errF <- learnErrors(opt$`filtF`, multithread=FALSE, nbases=opt$`nbases`)

   flog.info("plotting error model")
   p <- plotErrors(errF, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errF.png", sep="_")
   ggsave(filename, height=10, width=10, type="cairo")
   p <- plotErrors(errR, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errR.png", sep="_")
   ggsave(filename, type="cairo")

   # dereplication
   flog.info("dereplication")
   derepF <- derepFastq(opt$`filtF`, verbose=TRUE)
   derepR <- derepFastq(opt$`filtR`, verbose=TRUE)

   # sample inference
   flog.info("sample inference")
   dadaF <- dada(derepF,
                 err=errF,
		 multithread=FALSE,
		 OMEGA_A=opt$`omega-a`,
		 USE_QUALS=opt$`use-quals`,
		 USE_KMERS=opt$`use-kmers`,
		 KDIST_CUTOFF=opt$`kdist-cutoff`,
		 BAND_SIZE=opt$`band-size`,
		 GAP_PENALTY=opt$`gap-penalty`,
		 HOMOPOLYMER_GAP_PENALTY=opt$`homopolymer-ga-penalty`,
		 MIN_FOLD=opt$`min-fold`,
                 MIN_HAMMING=opt$`min-hamming`,
                 MIN_ABUNDANCE=opt$`min-abundance`,
		 MAX_CLUST=opt$`max-clust`,
		 MAX_CONSIST=opt$`max-consist`)

   dadaR <- dada(derepR,
                 err=errR,
		 multithread=FALSE,
		 OMEGA_A=opt$`omega-a`,
		 USE_QUALS=opt$`use-quals`,
		 USE_KMERS=opt$`use-kmers`,
		 KDIST_CUTOFF=opt$`kdist-cutoff`,
		 BAND_SIZE=opt$`band-size`,
		 GAP_PENALTY=opt$`gap-penalty`,
		 HOMOPOLYMER_GAP_PENALTY=opt$`homopolymer-ga-penalty`,
		 MIN_FOLD=opt$`min-fold`,
                 MIN_HAMMING=opt$`min-hamming`,
                 MIN_ABUNDANCE=opt$`min-abundance`,
		 MAX_CLUST=opt$`max-clust`,
		 MAX_CONSIST=opt$`max-consist`)

   # merge pairs - returns a dataframe
   flog.info("merging paired reads")
   mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE, propagateCol=c("n0", "n1", "nunq", "pval", "birth_type", "birth_pval", "birth_fold", "birth_ham", "birth_qave", "abundance"))


   # get diagnostic clustering data frame
   df.clustering <- data.frame(sequence=mergers$sequence,
                               F.abundance=mergers$F.abundance,
			       F.n0=mergers$F.n0,
                               F.n1=mergers$F.n1,
			       F.nunq=mergers$F.nunq,
			       F.pval=mergers$F.pval,
			       F.birth_pval=mergers$F.birth_pval,
			       F.birth_fold=mergers$F.birth_fold,
			       F.birth_ham=mergers$F.birth_ham,
			       F.birth_qave=mergers$F.birth_qave,
                               R.abundance=mergers$R.abundance,
                               R.n0=mergers$R.n0,
                               R.n1=mergers$R.n1,
			       R.nunq=mergers$R.nunq,
			       R.pval=mergers$R.pval,
			       R.birth_pval=mergers$R.birth_pval,
			       R.birth_fold=mergers$R.birth_fold,
			       R.birth_ham=mergers$R.birth_ham,
			       R.birth_qave=mergers$R.birth_qave)


   clustering.filename <- paste0(opt$`outdir`, "/", sample.name, "_clustering.tsv")
   write.table(df.clustering, clustering.filename, sep="\t", row.names=F, quote=F)

   mergers <- data.frame(sequence=mergers$sequence,
   	                 abundance=mergers$abundance,
			 forward=mergers$forward,
			 reverse=mergers$reverse,
			 nmatch=mergers$nmatch,
			 nmismatch=mergers$nmismatch,
			 nindel=mergers$nindel,
			 prefer=mergers$prefer,
			 accept=mergers$accept)

   flog.info("reoving sequences where accept==FALSE")
   # just keep accepted merged sequences
   mergers <- mergers[mergers$accept==TRUE,]

   # remove chimeras
   flog.info("removing chimeric sequences")
   mergers.nochim <- removeBimeraDenovo(mergers, method="consensus")

   # write out seq table
   flog.info("writing results")
   outfile <- paste(opt$`outdir`, sample.name, sep="/")
   outfile <- paste(outfile, "_seq_abundance.tsv", sep="")
   write.table(mergers.nochim, file=outfile, sep="\t", quote=F, row.names=F)
   
   # get summaries of reads passing each stage
   flog.info("writing summaries")
   getN <- function(x) sum(getUniques(x))
   track <- data.frame(cbind(getN(dadaF), getN(dadaR), getN(mergers), sum(mergers.nochim$abundance)))
   colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
   track$sample <- sample.name

   summary.outfile <- paste(opt$`outdir`, sample.name, sep="/")
   summary.outfile <- paste(summary.outfile, "_summary.tsv", sep="")
   write.table(track, summary.outfile, sep="\t", row.names=F, quote=F)
}