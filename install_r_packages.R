# Little script to check if R packages are installed and
# install if not

packages <- c("dplyr",  
              "ggplot2",                                                                     "gplots",                                                                      "gridExtra",                                                                   "gtools",                                                                      "rmarkdown",                                                                   "optparse",                                                                    "plotly",                                                                      "RColorBrewer",                                                                "RcppParallel",                                                                "reshape")


# This was taken from
# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
# excluding the loading part

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
  }
)