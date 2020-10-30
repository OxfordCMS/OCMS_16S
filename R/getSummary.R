# get summary statistics from a list of vectors

getSummary <- function(values = list(), func=mean){

    summary.list = list()
    for (i in 1:length(values)){
        summary.list[[i]] <- func(values[[i]])
    }
    summary.list <- unlist(summary.list)
    return(summary.list)
   }


