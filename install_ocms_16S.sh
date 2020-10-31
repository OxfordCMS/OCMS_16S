
# r dependencies
echo "installing R dependencies"
conda install -c r r-ggplot2 \
      r-plotly \
      r-gridextrar \
      r-dplyr \
      r-reshape \
      r-gplots \
      r-data.table \
      r-optparse \
      r-cppparallel

echo "installing dada2"
conda install -c bioconda bioconductor-dada2

echo "installing cgat-core"
conda install -c bioconda cgat-core

echo "installing OCMS 16S"
python setup.py install
