###################################################
###################################################
# Filtering and trimming options for dada2
###################################################
###################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dada2"))
suppressPackageStartupMessages(library("futile.logger"))

# make options list
option_list <- list(
               make_option(c("-d", "--directory"), default=".",
                           help="directory of fastq files [default %default]"),
               make_option(c("-f", "--filtered-directory"), default="filtered",
                           help="directory for filtered fastq files [default %default]"),
	       make_option(c("-p", "--paired"), action="store_true", default=FALSE,
                           help="is it paired-end data [default %default]"),
	       make_option(c("-n", "--maxN"), default=0,
                           help="maxN parameter [default %default]"),
               make_option(c("-e", "--maxEE"), default="2,2",
                           help="maximum number of expected errors [default %default]"),
               make_option(c("--truncQ"), default=2,
                           help="truncate reads at the first instance of a quality score less than or equal to truncQ [default %default]"),
               make_option(c("--truncLen"), default="250,250",
                           help="truncate  reads  after truncLen bases [default %default]"))

# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################
opt <- parse_args(OptionParser(option_list=option_list))

###########################
# helper functions
###########################

getInputFastq <- function(directory, paired){

	         if (paired){
		    fnFs <- sort(list.files(directory, pattern=".fastq.1", full.names = TRUE))
		    fnRs <- sort(list.files(directory, pattern=".fastq.2", full.names = TRUE))}
		 else{
                    fnFs <- sort(list.files(directory, pattern=".fastq.1", full.names = TRUE))
		    fnRs <- NA}

		 return(list(fnFs, fnRs))
		 }

###########################
###########################
###########################

getSampleNames <- function(fnFs){

	       sample.names <- gsub(".fastq.1", "", basename(fnFs))
	       return(sample.names)
	       }

###########################
###########################
###########################

splitArg <- function(arg){

	 # split input list argument
	 arg <- strsplit(arg, ",")
	 return(unlist(arg))
	 }
	 
###########################
###########################
###########################

# list options
truncLen <- as.numeric(splitArg(opt$`truncLen`))
maxEE <- as.numeric(splitArg(opt$`maxEE`))

# select input fastq files
input.files <- getInputFastq(opt$`directory`, opt$`paired`)

# input files
fnFs <- unlist(input.files[1])
fnRs <- unlist(input.files[2])

# sample names
sample.names <- getSampleNames(fnFs)

# filtering
filtFs <- file.path(opt$`filtered-directory`, paste0(sample.names, ".fastq.1.gz"))

if (!(is.na(fnRs))){
   filtRs <- file.path(opt$`filtered-directory`, paste0(sample.names, ".fastq.2.gz"))
   }

flog.info(paste0("filtering reads and writing to ", opt$`filtered-directory`))
if (!(is.na(fnRs))){
   out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        truncLen=truncLen,
                        maxN=opt$`maxN`,
			maxEE=maxEE, truncQ=opt$`truncQ`, rm.phix=TRUE,
	                compress=TRUE, multithread=TRUE)
   flog.info(paste0(paste0("writing summary to ", opt$`filtered-directory`), "/summary.tsv"))
   outfile <- paste(opt$`filtered-directory`, "summary.tsv", sep="/")
   write.table(out, file=outfile, sep="\t")
   }

if (is.na(fnRs)){
   out <- filterAndTrim(fnFs, filtFs,
                        truncLen=truncLen,
                        maxN=opt$`maxN`,
			maxEE=maxEE, truncQ=opt$`truncQ`, rm.phix=TRUE,
	                compress=TRUE, multithread=TRUE)
   flog.info(paste0(paste0("writing summary to ", opt$`filtered-directory`), "/summary.tsv"))
   outfile <- paste(opt$`filtered-directory`, "summary.tsv", sep="/")
   write.table(out, file=outfile, sep="\t")
   }

flog.info("DONE")




