=======================================================
run dada2 amplicon sequencing analysis
=======================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
--------

The purpose of this pipeline is to take amplicon sequencing data and run it through
the DADA2 pipeline - this includes filtering and trimming reads, merging and clustering
reads and assigning reads to taxa.

Dependencies
-------------

In order to run this pipeline you will need to install the following packages:

* CGAT tools. The documentation can be found `here`_  

.. _here: https://github.com/cgat-developers/cgat-flow

* NGSKit. The code can be found in this `github repository`_.

.. _github repository: https://github.com/nickilott/NGSKit

Installation
-------------

You should install CGAT tools using the installation script that is provided, this will
install everything you need into a conda environment. 

To activate the enviroment do the following::

    source </full/path/to/folder/without/trailing/slash>/conda-install/etc/profile.d/conda.sh
    conda activate base
    conda activate cgat-f

    # finally, please run the cgatflow command-line tool to check the installation:
    cgatflow --help

.. warning::

    Make sure that you are installing CGAT into a clean environment i.e. with no modules
    loaded from shared directories. You can do something like module purge to remove loaded
    modules.

However it does not come with dada2 installed. Dada2 requires RcppParallel to be installed so
first, so once you have activated the CGAT environment install this::

    conda install -c conda-forge r-rcppparallel 

Then you should install the dada2 R package as follows::

    conda install -c bioconda bioconductor-dada2

This should install dada2 so that when you start R you can simply type::

    library("dada2")

You will in all likelihood also have to install optparse for argument parsing in the R scripts::

    conda install -c conda-forge r-optparse

Once you have successfully installed CGAT, dada2 and optparse, you should then download the NGSKit
repository - this contains the dada2 pipeline (/pipelines/pipeline_dada2.py). In a development
folder run::
  
    git clone git@github.com:nickilott/NGSKit.git

It is possible that the modules that are required for the pipeline to run will have to be
placed into your PYTHONPATH::

    PYTHONPATH=$PYTHONPATH:<path-to-NGSKit>
    export $PYTHONPATH

Input files
------------

The input is a directory of fastq formatted files. These should be placed in the directory in
which you wish to run the pipeline. They must be of the format <name>.fastq.1.gz for single-end
data and two files <name>.fastq.1.gz and <name>.fastq.2.gz for paired-end data.

Parameterisation
------------------

There is a pipeline.yml file that specifies the parameters to be passed to dada2 and reporting
functions in the pipeline. This is located in <path-to-NGSKit>/pipelines/pipeline_dada2/pipeline.yml
and needs to be copied to your working directory. The parameters should be changed to suitable values.
and are explained in the dada2 `tutorial`_. 

.. _tutorial: https://benjjneb.github.io/dada2/tutorial.html

Output files
-------------

The main output file of the pipeline is the counts matrix that consists of amplicon sequence variants
and their abundance in each sample. The pipeline assigns taxonomy to each ASV and this is incorporated
into the ASV name in the resulting file. It is of the form:

+---------------------------------------------------------------------+---------+----------+
|test_id                                                              | Sample1 | Sample2  |
+---------------------------------------------------------------------+---------+----------+
|ASV1:p__phylum1;c__class1;o__order1;f__family1;g__genus1;s__species1 | 1000    | 1239     |
+---------------------------------------------------------------------+---------+----------+
|ASV2:p__phylum2;c__class2;o__order2;f__family2;g__genus2;s__species2 | 500     | 10       |
+---------------------------------------------------------------------+---------+----------+
|ASV3:p__phylum3;c__class3;o__order3;f__family3;g__genus3;s__species3 | 1000    | 2300     |
+---------------------------------------------------------------------+---------+----------+


The purpose of this output file is that it can be taken forward in a easy fashion to look at differential
abundance using software such as DESeq2 and this will be done on a per ASV level. If you wish to perform
analysis on counts that have been summed over taxa at a particular taxonomic level you can use the following
output files:

* taxonomy_abundances.dir/phylum_abundances.tsv
* taxonomy_abundances.dir/class_abundances.tsv
* taxonomy_abundances.dir/order_abundances.tsv
* taxonomy_abundances.dir/family_abundances.tsv
* taxonomy_abundances.dir/genus_abundances.tsv
* taxonomy_abundances.dir/species_abundances.tsv

Running the pipeline
---------------------

Running the pipeline should be as simple as::

    python <path-to-NGSKit>/pipelines/pipeline_dada2.py -v5 -p50 make full

Set p to the number of processes that you want to run - this corresponds to the number of the samples
that you want to run (unless it is a huge number and exceeds a reasonable number to run on the cluster).

Reporting
----------

Once the pipeline has finished you can build a report by running::

    python <path-to-NGSKit>/pipelines/pipeline_dada2.py -v5 -p50 make build_report
