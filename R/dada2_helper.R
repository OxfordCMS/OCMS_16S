###############################################
###############################################
###############################################
# Manipulating dada2 output
###############################################
###############################################
###############################################

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
  		           c <- gsub("ASV[0-9]*:p__.*;", "",  c)
                           d <- aggregate(dat, by=list(c), sum)}
                       else if (level == "order"){
                           o <- gsub(";f__.*", "", rownames(dat))
		           o <- gsub("ASV[0-9]*:p__.*;c__.*;", "",  o)
                           d <- aggregate(dat, by=list(o), sum)}
		       else if (level == "family"){
                           f <- gsub(";g__.*", "", rownames(dat))
		           f <- gsub("ASV[0-9]*:p__.*;c__.*;o__.*;", "",  f)
                           d <- aggregate(dat, by=list(f), sum)}
		       else if (level == "genus"){
		           g <- gsub(";s__.*", "", rownames(dat))
		           g <- gsub("ASV[0-9]*:p__.*;c__.*;o__.*;f__.*;", "",  g)
                           d <- aggregate(dat, by=list(g), sum)}
		       else if (level == "species"){
		           s <- gsub(".*;", "", rownames(dat))
		           g <- gsub(";s__.*", "", rownames(dat))
		           g <- gsub("ASV[0-9]*:p__.*;c__.*;o__.*;f__.*;", "",  g)
			   s <- gsub("g__", "s__", gsub("s__", "", paste(g, s, sep="_")))
                           d <- aggregate(dat, by=list(s), sum)}
                   }
		   else{
		   stop("level not one of phylum, class, order, family, genus, species")
		   }
		   colnames(d)[colnames(d) == "Group.1"] <- "Taxon"
		   return(d)
		   }

		   