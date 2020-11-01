
# r dependencies
echo "installing R dependencies"
conda install -c condo-forge python=3.7.1
conda install -c conda-forge r-base=3.6.1
conda install -c condo-forge r-rcppparallel
conda install -c condo-forge r-optparse
conda install -c conda-forge r-ggplot2 \
      r-plotly \
      r-gridextra \
      r-dplyr \
      r-reshape \
      r-gplots \
      r-data.table \

echo "installing bioconda packages"
conda install -c cgat-core bioconductor-dada2

echo "installing OCMS 16S"
python setup.py install
