###################################################
###################################################
# Error model, de-replication and sample inference
###################################################
###################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dada2"))
suppressPackageStartupMessages(library("futile.logger"))
suppressPackageStartupMessages(library("dplyr"))

# make options list
option_list <- list(
               make_option(c("--seqfile"),
                           help="file with sequence and abundance"),
               make_option(c("--training-set"),
                           help="training set for assigning taxonomy"),
               make_option(c("--species-file"),
                           help="file for assigning at the species level"),
               make_option(c("-o", "--outfile"),
                           help="outfile containing taxa and abundances")
			   )

# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################

opt <- parse_args(OptionParser(option_list=option_list))

# read sequence to abundance file
seqs <- read.csv(opt$`seqfile`, header=T, stringsAsFactors=F, sep="\t")

# assign taxonomy
flog.info("assigning taxonomy")
taxa <- assignTaxonomy(seqs, opt$`training-set`)

# add species
flog.info("adding species assignments")
taxa <- addSpecies(taxa, opt$`species-file`)
taxa <- as.data.frame(unlist(taxa))

# merge the data
flog.info("merging data")
all.data <- merge(taxa, seqs, by.x="row.names", by.y="sequence", all.x=T, all.y=T)
colnames(all.data) <- gsub("Row.names", "sequence", colnames(all.data))

# write data
flog.info("writing output")
write.table(all.data, file=opt$`outfile`, sep="\t", quote=F, row.names=F)

