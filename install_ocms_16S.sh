
# r dependencies
echo "installing R dependencies"
conda install -c r r-ggplot2 \
      r-plotly \
      r-gridextra \
      r-dplyr \
      r-reshape \
      r-gplots \
      r-data.table \
      r-optparse \

echo "installing bioconda packages"
conda install -c bioconda r-rcppparallel bioconductor-dada2 cgat-core

echo "installing OCMS 16S"
python setup.py install
