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
               make_option(c("-i", "--infile"), default=NA, type="character",
                           help="input fastq file [default %default]"),
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
                           help="truncate  reads  after truncLen bases [default %default]"),
               make_option(c("--trimLeft"), default="0,0",
                           help="trim left sequence (primers) [default %default]"))


# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################
opt <- parse_args(OptionParser(option_list=option_list))

###########################
# helper functions
###########################

getInputFastq <- function(infile, paired){

	         if (paired){
		    fnF <- infile
		    fnR <- gsub(".fastq.1", ".fastq.2", fnF)}
		 else{
                    fnF <- infile
		    fnR <- NA}

		 return(list(fnF, fnR))
		 }

###########################
###########################
###########################

getSampleName <- function(fnF){

	       sample.name <- gsub(".fastq.1", "", basename(fnF))
	       return(sample.name)
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
trimLeft <- as.numeric(splitArg(opt$`trimLeft`))

# select input fastq files

input.files <- getInputFastq(opt$`infile`, opt$`paired`)

# input files
fnF <- unlist(input.files[1])
fnR <- unlist(input.files[2])

# sample names
sample.name <- getSampleName(fnF)

# filtering
filtF <- file.path(opt$`filtered-directory`, paste0(sample.name, ".fastq.1.gz"))

if (!(is.na(fnR))){
   filtR <- file.path(opt$`filtered-directory`, paste0(sample.name, ".fastq.2.gz"))
   }

flog.info(paste0("filtering reads and writing to ", opt$`filtered-directory`))
if (!(is.na(fnR))){

   out <- filterAndTrim(fnF, filtF, fnR, filtR,
                        truncLen=truncLen,
			trimLeft=trimLeft,
                        maxN=opt$`maxN`,
			maxEE=maxEE, truncQ=opt$`truncQ`, rm.phix=TRUE,
	                compress=TRUE, multithread=FALSE)
   flog.info(paste0(paste0("writing summary to ", opt$`filtered-directory`), "/summary.tsv"))
   out <- as.data.frame(out)
   out$sample <- sample.name
   outfile <- paste(opt$`filtered-directory`, sample.name, sep="/")
   outfile <- paste(outfile, "summary.tsv", sep="_")
   write.table(out, file=outfile, sep="\t", row.names=F, quote=F)
   }

if (is.na(fnR)){
   out <- filterAndTrim(fnF, filtF,
                        truncLen=truncLen,
			trimLeft=trimLeft,
                        maxN=opt$`maxN`,
			maxEE=maxEE, truncQ=opt$`truncQ`, rm.phix=TRUE,
	                compress=TRUE, multithread=FALSE)
   flog.info(paste0(paste0("writing summary to ", opt$`filtered-directory`), "/summary.tsv"))
   out <- as.data.frame(out)
   out$sample <- sample.name
   outfile <- paste(opt$`filtered-directory`, sample.name, sep="/")
   outfile <- paste(outfile, "summary.tsv", sep="_")
   write.table(out, file=outfile, sep="\t", row.names=F, quote=F)
   }

flog.info("DONE")




