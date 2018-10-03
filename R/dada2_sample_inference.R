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
	       make_option(c("-n", "--nreads"), default=1000000,
                           help="number of reads to learn error model [default %default]"),
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

flog.info(gsub("nreads", as.character(opt$`nreads`), "learning errors on the first nreads reads"))

if (is.na(opt$`filtR`)){

   # error model inference
   errF <- learnErrors(opt$`filtF`, multithread=TRUE, nreads=opt$`nreads`)
   flog.info("plotting error model")
   p <- plotErrors(errF, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errF.pdf", sep="_")
   ggsave(filename, height=10, width=10)

   # dereplication
   flog.info("dereplication")
   derepF <- derepFastq(opt$`filtF`, verbose=TRUE)

   # sample inference
   flog.info("sample inference")	
   dadaF <- dada(derepF, err=errF, multithread=TRUE)

   # return a dataframe - I find this more intuitive
   dadaF.df <- as.data.frame(dadaF$denoised)
   dadaF.df$sequence <- rownames(dadaF.df)
   dadaF.df <- data.frame(sequence=dadaF.df$sequence, abundance=dadaF.df[,1])

   # remove chimeras
   flog.info("removing chimeric sequences")
   dadaF.df <- removeBimeraDenovo(dadaF.df, method="consensus")

   # write out table
   flog.info("writing results")
   outfile <- paste(opt$`outdir`, sample.name, sep="/")
   outfile <- paste(outfile, "_seq_abundance.tsv.gz", sep="")
   write.table(dadaF.df, file=outfile, sep="\t", quote=F, row.names=F)

} else {

   # error model learning
   errR <- learnErrors(opt$`filtR`, multithread=TRUE, nreads=opt$`nreads`)
   errF <- learnErrors(opt$`filtF`, multithread=TRUE, nreads=opt$`nreads`)

   flog.info("plotting error model")
   p <- plotErrors(errF, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errF.pdf", sep="_")
   ggsave(filename, height=10, width=10)
   p <- plotErrors(errR, nominalQ=TRUE)
   filename <- paste(paste(directory, sample.name, sep="/"), "errR.pdf", sep="_")
   ggsave(filename, height=10, width=10)

   # dereplication
   flog.info("dereplication")
   derepF <- derepFastq(opt$`filtF`, verbose=TRUE)
   derepR <- derepFastq(opt$`filtR`, verbose=TRUE)

   # sample inference
   flog.info("sample inference")
   dadaF <- dada(derepF, err=errF, multithread=TRUE)
   dadaR <- dada(derepR, err=errR, multithread=TRUE)

   # merge pairs
   mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE)

   # write out seq table
   flog.info("writing results")
   outfile <- paste(opt$`outdir`, sample.name, sep="/")
   outfile <- paste(outfile, "_seq_abundance.tsv.gz", sep="")
   write.table(mergers, file=outfile, sep="\t", quote=F, row.names=F)
   }
   
   