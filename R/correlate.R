
library(dplyr)

correlateMatrixWithVector <- function(mat, vector, method="spearman"){
  
  cors <- list()
  for (i in 1:nrow(mat)){
    res <- cor.test(unlist(mat[i,]), vector, method=method)
    res <- data.frame(cor=res$estimate, p.value=res$p.value)
    cors[[i]] <- res
  }
  result <- bind_rows(cors)
  rownames(result) <- rownames(mat)
  result$padj <- p.adjust(result$p.value)
  return(result)
}




