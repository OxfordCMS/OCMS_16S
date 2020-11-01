suppressPackageStartupMessages(library("optparse"))
library(rmarkdown)

option_list <- list(
               make_option(c("-i", "--infile"), default=NA, type="character",
                           help="input fastq file [default %default]"))

# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################
opt <- parse_args(OptionParser(option_list=option_list))

render(opt$`infile`, output_format="html_document")