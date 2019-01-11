
rowScale <- function(mat){

	 mat.s <- data.frame(t(apply(mat, 1, scale)))
	 return(mat.s)
	 }

runKmeans <- function(mat, k){

          set.seed(1984)
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
