# OCMS_16S - amplicon sequencing analysis

## Purpose

This is a pipeline that is built using the
[cgat-core](https://github.com/cgat-developers/cgat-core) framework. The
purpose of the pipeline is to run dada2 processing of amplicon
sequencing data either on a compute cluster or locally. The pipeline
consists of a number of wrapper scripts in R that are executed on the
commandline and as such there is no requirement for R coding by the
user. The hope is that the pipeline provides an accessible and
user-friendly interface to produce reproducble results from an amplicon
sequencing study.

A useful feature of this pipeline is that it outputs an sqlite database
that serves as input to
[OCMSlooksy](https://github.com/OxfordCMS/OCMSlooksy), an R/Shiny
application designed for downstream analysis and visualisation of
amplicon sequencing data.

You should familiarise yourself with the
[dada2](https://benjjneb.github.io/dada2/tutorial.html) workflow before
running the pipeline to ensure that you understand the parameterisation.

## Installation

The pipeline depends on having a number of python and R libraries
installed. There are various ways that you can install and use the
pipeline and these are outlined below.

### conda

It is possible to set up a conda environment that is compatible with our
dada2 pipeline that has all of the relevant dependencies installed. The
following steps outline how to install the pipeline using conda to
manage the dependencies.

1.  Download and install
    [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
    if you don\'t already have it.

2.  Clone the OCMS_16S repository:

        git clone git@github.com:OxfordCMS/OCMS_16S.git

3.  Change into OCMS_16S directory:

        cd OCMS_16S

4.  Create ocms_16s conda environment:

        conda env create -f envs/environment_ocms_16s.yaml

This can often be slow as conda tries to solve environment. If this
happens you should try using mamba instead:

    conda install -c conda-forge mamba

    mamba env create -f envs/environment_ocms_16s.yaml

5.  Activate the newly created environment:

        conda activate ocms_16s

6.  Install ocms_16s dada2:

        python setup.py install

### Python virtual environment

You can use a python virtual environment to install python dependencies
whilst relying on system-wide software for non-python dependencies (R
etc.).

1.  Create a python virtual environment with your chosen python install
    (call the environment what you like. We call it ocms_pipeline as an
    example):

        python -m venv ocms_pipeline

2.  Activate virtualenv:

        source ocms_pipeline/bin/activate

3.  Install cgat-core:

        # first install some dependencies
        pip install apsw gevent numpy==1.19.5 pandas paramiko pep8 pytest pytest-pep8 drmaa pyyaml ruffus setuptools six sqlalchemy

        pip install cgatcore

4.  Make sure you have other R dependencies installed:

        dplyr
        ggplot2
        gplots
        gridextra
        gtools
        rmarkdown
        optparse
        plotly
        rcolorbrewer
        rcppparallel
        reshape

You can install these one-by-one or using the script,
install_r\_packages.R, that is present in the top-level directory of the
OCMS_16S repository. Start R and type:

    source("install_r_packages.R")

5.  Install dada2 as described
    [here](https://benjjneb.github.io/dada2/dada-installation.html)

6.  Install OCMS_16S:

        pip install ocms-16s

## Running the pipeline

The pipeline runs using fastq files as input. It also requires
parameters to be set that will be passed to the dada2 R scripts.

### Input files

The input is a directory of fastq files. These should be placed in the
directory in which you wish to run the pipeline. They must be of the
format \<name\>.fastq.1.gz for single-end data and two files
\<name\>.fastq.1.gz and \<name\>.fastq.2.gz for paired-end data.

For the pipeline to run succesfully you will also need to have
downloaded relevant [dada2
databases](https://benjjneb.github.io/dada2/training.html) and point to
them in the pipeline.yml parameters file as described in the next
section.

### Parameterisation

The parameters for dada2 processing are specified in the pipeline.yml
file. To create this file, move into the directory containing the fastq
files that you wish to process and type:

    ocms_16s dada2 config

This will create the pipeline.yml file in the current working directory
which you can edit using your favourite text editor. The parameters are
provided in a standard yaml format as outlined below:

    # specify whether data are paired or single end. The
    # pipeline will pick up whether this is true but being
    # explicit here is helpful
    paired: 1

    # dada2 parameters
    trim:

        # parameters used for trimming reads. If the data are
        # paired-end then you need to specify 2 values for
        # maxee, truncLen and trimLeft. These parameters must be specified
        maxn: 0
        maxee: 2,2
        truncq: 2
        trunclen: 250,160
        trimleft: 0,0

    sample_inference:

        # parameters for sample inference. This includes
        # error learning, de-replication, merging (if paired) and
        # sample inference.

        # number of reads to use (per sample) to estimate error
        # model
        nbases: 10000000

        # additional options
        options: ''

    taxonomy:

        memory: 10G

        # assigning taxonomy
        taxonomy_file: RefSeq-RDP16S_v2_May2018.fa.gz

        # This is the file that is used for the addSpecies function in
        # dada2 for exact matching and species assignment. It must therefore
        # be derived from the same database used as taxonomy_file above
        species_file: silva_species_assignment_v132.fa.gz

    report:
        # whether to run diagnostics report. This is only necessary if after the
        # main report is built you want to get into more regarding the specifics of
        # how dada2 processed sequences. Specify as 1 if you wish to run it
        diagnostics:

    database:
        # name of the output database. This is a database that is built to
        # be compatible with the OCMSlooksy.
        name: output_db

The majority of the parameters correspond to the dada2 arguments to the
various functions in the dada2 package.

### Getting help on pipeline tasks

The pipeline is run using a simple commandline interface. You can view
the tasks that are going to be run by using the \'show\' command. In the
directory that you plan to run the pipeline:

    ocms_16s dada2 show full

This will print out the tasks that are going to be run:

    ----------------------------------------------------
    Tasks which will be run:

    Task = "mkdir('tree.dir')   before pipeline_dada2.buildTree "
    Task = "mkdir('abundance.dir')   before pipeline_dada2.runSampleInference "
    Task = "mkdir('filtered.dir')   before pipeline_dada2.filterAndTrim "
    Task = 'pipeline_dada2.filterAndTrim'
    Task = 'pipeline_dada2.runSampleInference'
    Task = 'pipeline_dada2.mergeAbundanceTables'
    Task = "mkdir('taxonomy.dir')   before pipeline_dada2.assignTaxonomy "
    Task = 'pipeline_dada2.assignTaxonomy'
    Task = 'pipeline_dada2.addUniqueIdentifiers'
    Task = 'pipeline_dada2.mergeTaxonomyTables'
    Task = 'pipeline_dada2.buildDefinitiveTable'
    Task = 'pipeline_dada2.buildTree'
    Task = "mkdir('taxonomy_abundances.dir')   before pipeline_dada2.splitTableByTaxonomicLevels "
    Task = 'pipeline_dada2.splitTableByTaxonomicLevels'
    Task = 'pipeline_dada2.full'
    ________________________________________
    # 2021-11-25 21:52:46,850 INFO job finished in 0 seconds at Thu Nov 25 21:52:46 2021 --  1.59  1.52  0.00  0.02 -- 1cae61fa-de0c-4b85-86b0-38dfd964c155

There are often numerous parameters that can be passed to dada2
functions. The most commone parameters that need to be changed are
explicitly stated in the pipeline.yml. However additional options can be
specified and these are commadline options to the various R scripts. You
can view these parameters by running the \'help\' script. For example:

    ocms_16s help --sampleInference

This will provide the possible options that can be passed to the
runSampleInference task via the pipeline.yml.

Once you have set the parameters, the pipeline should be simple to run.
You can run the pipeline locally or on a compute cluster in order to
maximise parallelisation that is afforded by using cgat-core workflow
management.

### Running the pipeline locally

In general the pipeline can be run using the following command:

    ocms_16s dada2 make full -v5 -p100

where -v specifies the verbosity level of the logging output and -p
specifies the number of processes you want to lauch per task e.g if you
want to process 100 samples then specifiy -p100 and each sample will be
processed in parallel and data combined in the final output tables.

If you want to run it locally on a laptop you will need access to a
unix-like operating system (e.g. Mac). You must specify the \--local
flag:

    ocms_16s dada2 make full -v5 -p1 --local

specifying -p as the number of processors you have available.

### Running the pipeline on a cluster

The best way to maximise the utility of the pipeline is to run it on a
high performance cluster - allowing you to parallelise sample
processing.To run on a cluster you will have to have a .cgat.yml file in
your home directory that specifies the queue manager, queue to use etc.
An example is belo:

    cluster:
        queue_manager: <slurm|sge|pbstorque>
        parallel_environment: <pe name>
        queue: <queue_name>

You will also need to make sure that the pipeline has access to the
drmaa library so it\'s best to set this as an environmental variable in
your \~/.bashrc:

    export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so

Once set up you should be able to run:

    ocms_16s dada2 make full -v5 -p100

As the pipeline runs, logging information will be printed to the screen
and also saved in the file pipeline.log. This file is useful to inspect
if the pipeline crashes and you need to debug.

### Building a report

Once the pipeline has finished, there is opportunity to assess the dada2
processing results in an html report by running:

    ocms_16s dada2 make build_report

This will build the report, report.dir/report.html which you can
inspect.

### Transition to OCMSlooksy

OCMSlooksy is an R/Shiny application that enables users to inspect data
from this dada2 processing pipeline as well as perform statistical
analysis and visualisation. By running:

    ocms_16s dada2 make build_db

you will build an sqlite database that contains all of the outputs
neccessary to load into OCMSlooksy. The database will be named according
to the specification in the pipipeline.yml. In the example above it
would be called \'output_db\' and this would be present in the current
working directory.

### Other output files

OCMS_16S will also output flat files that can be used for downstream
analysis. The main output file of the pipeline is the counts matrix that
consists of amplicon sequence variants and their abundance in each
sample. The pipeline assigns taxonomy to each ASV and this is
incorporated into the ASV name in the resulting file. It is of the form:

  ---------------------------------------------------------------------------------- --------- ---------
  test_id                                                                            Sample1   Sample2

  ASV1:p\_\_phylum1;c\_\_class1;o\_\_order1;f\_\_family1;g\_\_genus1;s\_\_species1   1000      1239

  ASV2:p\_\_phylum2;c\_\_class2;o\_\_order2;f\_\_family2;g\_\_genus2;s\_\_species2   500       10

  ASV3:p\_\_phylum3;c\_\_class3;o\_\_order3;f\_\_family3;g\_\_genus3;s\_\_species3   1000      2300
  ---------------------------------------------------------------------------------- --------- ---------

This file is created as abundance.dir/taxa_abundances.tsv.

The purpose of this output file is that it can be taken forward in a
easy fashion to look at differential abundance using software such as
DESeq2 and this will be done on a per ASV level. If you wish to perform
analysis on counts that have been summed over taxa at a particular
taxonomic level you can use the following output files:

-   taxonomy_abundances.dir/phylum_abundances.tsv
-   taxonomy_abundances.dir/class_abundances.tsv
-   taxonomy_abundances.dir/order_abundances.tsv
-   taxonomy_abundances.dir/family_abundances.tsv
-   taxonomy_abundances.dir/genus_abundances.tsv
-   taxonomy_abundances.dir/species_abundances.tsv

# Acknowledgements

This pipeline is based off of a lot of work that has gone before it. It
is basically a wrapper for dada2 functionality and so if you use the
pipeline in a publication please remember to cite the dada2
[paper](https://www.nature.com/articles/nmeth.3869). The [cgat-core
framework](https://f1000research.com/articles/8-377x) is of course
another important tool that has enabled the development of this
pipeline.
