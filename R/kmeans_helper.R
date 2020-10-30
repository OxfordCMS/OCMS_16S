
rowScale <- function(mat){

	 mat.s <- data.frame(t(apply(mat, 1, scale)))
	 return(mat.s)
	 }

runKmeans <- function(mat, k){

  	kmeans.obj <- kmeans(mat, k)
	  return(kmeans.obj)
	  }

getWithinClusterVar <- function(kmeans.obj){

		  return(kmeans.obj$tot.withinss)
		  }

getClusters <- function(kmeans.obj){

	    return(data.frame(kmeans.obj$cluster))
	    }

elbow <- function(mat, ks=seq(1, 10, 1)){

      vars <- c()
      for (i in ks){
      	  kmeans.obj <- runKmeans(mat, i)
          var <- getWithinClusterVar(kmeans.obj)
	    vars <- append(vars, var)
      }
      df <- data.frame(k=ks, within.var=vars)
      return(df)
      }

getClusterGenes <- function(clusters, cluster=1){

		c.genes <- rownames(clusters)[clusters$kmeans.obj.cluster == cluster]
		return(c.genes)
		}

runSilhouette <- function(mat, ks=c(1:10), clustering_method="manhattan"){
  
  ave_sil_widths <- c()
  for (i in ks){
    if (i == 1){
      ave_sil_width=0
      }
    else{
    kmeans.obj <- runKmeans(mat, i)
    clustering <- kmeans.obj$cluster
    sil <- silhouette(clustering, dist(mat, method=clustering_method))
    ave_sil_width <- mean(sil[,"sil_width"])
      }
    ave_sil_widths <- append(ave_sil_widths, ave_sil_width)
  }
  df <- data.frame(k=ks, ave_sil_width=ave_sil_widths)
  return(df)
  }


