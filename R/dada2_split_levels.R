###################################################
###################################################
# split output file at each taxonomic level
###################################################
###################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("futile.logger"))

# make options list
option_list <- list(
               make_option(c("-i", "--infile"), default=NA, type="character",
                           help="input data table from dada2 pipeline [default %default]"),
               make_option(c("--outdir"), default="abundance.dir",
                           help="output for files [default %default]"))

# suppress warning messages
options(warn=-1)

###########################
# get command line options
###########################
opt <- parse_args(OptionParser(option_list=option_list))

buildPerLevelTable <- function(taxa_abundances, level="phylum"){

		   possible.levels = c("phylum",
		   		       "class",
				       "order",
				       "family",
				       "genus",
				       "species")
			
		   if (level %in% possible.levels){
		       if (level == "phylum"){
           		   p <- gsub(";c__.*", "", rownames(dat))
		           p <- gsub("ASV[0-9]*:", "", p)
   		           d <- aggregate(dat, by=list(p), sum)}
                       else if (level == "class"){
			   c <- gsub(";o__.*", "", rownames(dat))
  		           c <- gsub("ASV[0-9]*:", "",  c)
                           d <- aggregate(dat, by=list(c), sum)}
                       else if (level == "order"){
                           o <- gsub(";f__.*", "", rownames(dat))
		           o <- gsub("ASV[0-9]*:", "",  o)
                           d <- aggregate(dat, by=list(o), sum)}
		       else if (level == "family"){
                           f <- gsub(";g__.*", "", rownames(dat))
		           f <- gsub("ASV[0-9]*:", "",  f)
                           d <- aggregate(dat, by=list(f), sum)}
		       else if (level == "genus"){
		           g <- gsub(";s__.*", "", rownames(dat))
		           g <- gsub("ASV[0-9]*:", "",  g)
                           d <- aggregate(dat, by=list(g), sum)}
		       else if (level == "species"){
		           s <- gsub(".*;", "", rownames(dat))
		           g <- gsub(";s__.*", "", rownames(dat))
		           g <- gsub("ASV[0-9]*:", "",  g)
			   s <- paste(g, s, sep=";")
                           d <- aggregate(dat, by=list(s), sum)}
                   }
		   else{
		   stop("level not one of phylum, class, order, family, genus, species")
		   }
		   colnames(d)[colnames(d) == "Group.1"] <- "Taxon"
		   return(d)
		   }

dat <- read.csv(opt$`infile`, header=T, stringsAsFactors=F, sep="\t", row.names=1, quote="")

flog.info("writing per level tables")
p <- buildPerLevelTable(dat, level="phylum")
outfile <- paste(opt$`outdir`, "phylum_abundance.tsv", sep="/")
write.table(p, file=outfile, sep="\t", quote=F, row.names=F)

c <- buildPerLevelTable(dat, level="class")
outfile <- paste(opt$`outdir`, "class_abundance.tsv", sep="/")
write.table(c, file=outfile, sep="\t", quote=F, row.names=F)

o <- buildPerLevelTable(dat, level="order")
outfile <- paste(opt$`outdir`, "order_abundance.tsv", sep="/")
write.table(o, file=outfile, sep="\t", quote=F, row.names=F)

f <- buildPerLevelTable(dat, level="family")
outfile <- paste(opt$`outdir`, "family_abundance.tsv", sep="/")
write.table(f, file=outfile, sep="\t", quote=F, row.names=F)

g <- buildPerLevelTable(dat, level="genus")
outfile <- paste(opt$`outdir`, "genus_abundance.tsv", sep="/")
write.table(g, file=outfile, sep="\t", quote=F, row.names=F)

s <- buildPerLevelTable(dat, level="species")
outfile <- paste(opt$`outdir`, "species_abundance.tsv", sep="/")
write.table(s, file=outfile, sep="\t", quote=F, row.names=F)

