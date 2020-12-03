'''
=======================================================
ocms_16s dada2 - amplicon sequencing analysis 
=======================================================

:Author: Nick Ilott, Sandi Yen
:Release: $Id$
:Date: |today|
:Tags: Python

--------
Purpose
--------

This is a pipeline that is built using the `cgat-core`_ framework. The purpose of the pipeline is to run dada2 processing of amplicon sequencing data either on a compute cluster or locally. The pipeline consists of a number of wrapper scripts in R that are executed on the commandline and as such there is no requirement for R coding by the user. The hope is that the pipeline provides an accessible and user-friendly interface to produce reproducble results from an amplicon sequencing study.

You should familiarise yourself with the `dada2`_ workflow before running the pipeline to ensure that you understand the parameterisation.


.. _cgat-core: https://github.com/cgat-developers/cgat-core

.. _dada2: https://benjjneb.github.io/dada2/tutorial.html 

-------------
Installation
-------------

The pipeline depends on having a number of python and R libraries installed. The easiest way at the moment to ensure that you have an environment that is compatible with our dada2 pipeline is to create a conda environment that has all of the relevant dependencies installed. The following steps outline how to install our dada2 pipeline.

1. Download and install `conda`_ if you don't already have it.

.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

2. Clone the OCMS_16S repository::

    git clone git@github.com:OxfordCMS/OCMS_16S.git

3. Change into OCMS_16S directory::

    cd OCMS_16S

4. Create ocms_16s conda environment::

    conda env create -f envs/environment_ocms_16s.yaml

5. Activate the newly created environment::

    conda activate ocms_16s

6. Install ocms_16s dada2::

    python setup.py install

Now you are ready to go!

---------------------
Running the pipeline
---------------------

The pipeline runs using fastq files as input. It also requires parameters to be set that will be passed to the dada2 R scripts.

Input files
------------

The input is a directory of fastq formatted files. These should be placed in the directory in which you wish to run the pipeline. They must be of the format <name>.fastq.1.gz for single-end data and two files <name>.fastq.1.gz and <name>.fastq.2.gz for paired-end data.

For the pipeline to run succesfully you will also need to have downloaded relevant `dada2 databases`_ and point to them in the pipeline.yml parameters file as described in the next section.


.. _dada2 databases: https://benjjneb.github.io/dada2/training.html

Parameterisation
------------------

The parameters for dada2 processing are specified in the pipeline.yml file. To create this file, move into the directory containing the fastq files that you wish to process and type::

    ocms_16s dada2 config

This will create the pipeline.yml file in the current working directory which you can edit using your favourite text editor. The parameters are provided in a standard yaml format as outlined below::

    # specify whether data are paired or single end. The
    # pipeline will pick up whether this is true but being
    # explicit here is helpful
    paired:

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
        taxonomy_file: /gfs/mirror/dada2/RefSeq-RDP16S_v2_May2018.fa.gz

        # This is the file that is used for the addSpecies function in
        # dada2 for exact matching and species assignment. It must therefore
        # be derived from the same database used as taxonomy_file above
        species_file: /gfs/mirror/dada2/silva_species_assignment_v132.fa.gz

    report:
        # whether to run diagnostics report. This is only necessary if after the
        # main report is built you want to get into more regarding the specifics of
        # how dada2 processed sequences. Specify as 1 if you wish to run it
        diagnostics:

        # author and name of the project for reporting purposes
        author: Nick Ilott
        title: Title

    database:
        # name of the output database. This is a database that is built to
        # be compatible with the OCMSExplorer.
        name: output_db


The majority of the parameters correspond to the dada2 arguments to the various functions in the dada2 package. You should familiarise yourself with these.


Executing the pipeline
-----------------------

Once you have set the parameters, the pipeline should be simple to run. If you are running on the cluster you can type::

    ocms_16s dada2 make full -v5 -p100

where -v specifies the verbosity level of the logging output and -p specifies the number of processes you want to lauch per task e.g if you want to process 100 samples then specifiy -p100 and each sample will be processed in parallel and data combined in the final output tables. This will run through the dada2 workflow and produce the output files described in the next section. 

We hope that the pipeline is not restricted to those that do not have access to a cluster. Nevertheless, to run the pipeline on a laptop you will need access to a unix-like operating system (e.g. Mac). to run locally you can add the --local flag to the command::

    ocms_16s dada2 make full -v5 -p8 --local

specifying -p as the number of available processors you have on your machine.

As the pipeline runs, logging information will be printed to the screen and also saved in the file pipeline.log. This file is useful to inspect if the pipeline crashes and you need to debug.

Output files
-------------

The main output file of the pipeline is the counts matrix that consists of amplicon sequence variants and their abundance in each sample. The pipeline assigns taxonomy to each ASV and this is incorporated into the ASV name in the resulting file. It is of the form:

+---------------------------------------------------------------------+---------+----------+
|test_id                                                              | Sample1 | Sample2  |
+---------------------------------------------------------------------+---------+----------+
|ASV1:p__phylum1;c__class1;o__order1;f__family1;g__genus1;s__species1 | 1000    | 1239     |
+---------------------------------------------------------------------+---------+----------+
|ASV2:p__phylum2;c__class2;o__order2;f__family2;g__genus2;s__species2 | 500     | 10       |
+---------------------------------------------------------------------+---------+----------+
|ASV3:p__phylum3;c__class3;o__order3;f__family3;g__genus3;s__species3 | 1000    | 2300     |
+---------------------------------------------------------------------+---------+----------+

This file is created as abundance.dir/taxa_abundances.tsv.


The purpose of this output file is that it can be taken forward in a easy fashion to look at differential abundance using software such as DESeq2 and this will be done on a per ASV level. If you wish to perform analysis on counts that have been summed over taxa at a particular taxonomic level you can use the following output files:

* taxonomy_abundances.dir/phylum_abundances.tsv
* taxonomy_abundances.dir/class_abundances.tsv
* taxonomy_abundances.dir/order_abundances.tsv
* taxonomy_abundances.dir/family_abundances.tsv
* taxonomy_abundances.dir/genus_abundances.tsv
* taxonomy_abundances.dir/species_abundances.tsv


Reporting
----------

The pipeline also has a standard report that can be built using::

    ocms_16s dada2 build_report -v5

This will build the html report that can be found in report.dir/report.html and provides various pieces of information regarding the processing of the data through the dada2 workflow including number of reads kept during each procesing step as well as some basic taxonomy informationl.


---------------------
Downstream analysis
---------------------

We have an example of downstream analysis that can be performed in R that can be found `here`_.

.. _here: https://oxfordcms.github.io/OCMS-blog/bioinformatics/Example-16S-rRNA-Analysis



Code
-----


'''

# load modules
from ruffus import *
import cgatcore.experiment as E
import cgatcore.pipeline as P
import shutil
import logging as L
import os
import sys
import glob
import re
import ocms16S.PipelineDada2 as PipelineDada2

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# Initialize the pipeline
P.initialize()
PARAMS = P.get_params()

# scripts directory - R scripts for dada2 functions
PARAMS["rscriptsdir"] = os.path.join(os.path.dirname(PipelineDada2.__file__), "R")

###################################################
###################################################
###################################################

INFILES = glob.glob("*.fastq.1.gz")
OUTFILES = ["filtered.dir/" + x for x in INFILES]

# check the existence of paired files
if PARAMS["paired"] == 1:
    for infile in INFILES:
        second_in_pair = P.snip(infile, ".fastq.1.gz") + ".fastq.2.gz"
        assert os.path.exists(second_in_pair), "no second read for file {infile}: should be {second_in_pair}"
        trunclen = PARAMS["trim_trunclen"].split(",")
        maxee = PARAMS["trim_maxee"].split(",")
        assert len(trunclen) == 2, "specify 2 values only for paired data (truncLen)"
        assert len(maxee) == 2, "specify 2 values only for paired data (maxee)"

###################################################
###################################################
###################################################
        
@follows(mkdir("filtered.dir"))
@transform(INFILES, regex("(\S+).fastq.1.gz"), r"filtered.dir/\1.fastq.1.gz")
def filterAndTrim(infile, outfile):
    '''
    filter and trim reads based on user-defined
    parameters
    '''
    maxn=PARAMS["trim_maxn"]
    maxee=PARAMS["trim_maxee"]
    truncQ=PARAMS["trim_truncq"]
    truncLen=PARAMS["trim_trunclen"]
    trimLeft=PARAMS["trim_trimleft"]

    # make sure parameters are present
    assert maxn != "" and maxee != "" and truncQ != "" and truncLen != "" and trimLeft !="", \
    "must specify all parameters to filterAndTrim"
    
    if PARAMS["paired"] == 1:
        paired = "--paired"
        # unzip one by one
        infile_read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")
        infiles = [infile, infile_read2]
    else:
        paired = ""
        infiles = [infile]

    tmpdir = P.get_temp_dir()
    for inf in infiles:
        outtmp = os.path.join(tmpdir, inf.replace(".gz", ""))
        statement = '''gunzip -c %(inf)s > %(outtmp)s'''
        P.run(statement)

    # hackabout
    outtmp = os.path.join(tmpdir, [x for x in infiles if x.endswith(".fastq.1.gz")][0].replace(".gz", ""))
    statement = '''Rscript %(rscriptsdir)s/dada2_filter_and_trim.R
                           --infile=%(outtmp)s
                           %(paired)s
                           --maxN=%(maxn)s
                           --maxEE=%(maxee)s
                           --truncQ=%(truncQ)s
                           --truncLen=%(truncLen)s
                           --trimLeft=%(trimLeft)s
                           --filtered-directory=filtered.dir'''
    P.run(statement)
    summary = open(outfile.replace(".fastq.1.gz",  "_summary.tsv"))
    summary.readline()
    summary = summary.readline().split("\t")
    assert int(summary[1]) != 0, f"Filtered file {outfile} has no reads - try re-running with different filtering parameters"
        
    shutil.rmtree(tmpdir)
    
###################################################
###################################################
###################################################

@follows(mkdir("abundance.dir"))
@transform(filterAndTrim, regex("(\S+)/(\S+).fastq.1.gz"), r"abundance.dir/\2_seq_abundance.tsv")
def runSampleInference(infile, outfile):
    '''
    run dada2 sample inference. It goes through the following
    stages:
    1) learning error model
    2) de-replication
    3) merging reads (if paired)
    4) removing chimeras
    '''
    nbases = PARAMS["sample_inference_nbases"]
    options = PARAMS["sample_inference_options"]
    outdir = os.path.dirname(outfile)
    job_memory=PARAMS["sample_inference_memory"]

    if PARAMS["paired"] == 1:
        # get second read in pair
        infile_read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")

        # generate statement to execute
        statement = '''Rscript %(rscriptsdir)s/dada2_sample_inference.R
                       --filtF=%(infile)s
                       --filtR=%(infile_read2)s
                       --nbases=%(nbases)s
                       %(options)s
                       --outdir=%(outdir)s'''
    else:
        statement = '''Rscript %(rscriptsdir)s/dada2_sample_inference.R
                       --filtF=%(infile)s
                       --nbases=%(nbases)s
                       %(options)s
                       --outdir=%(outdir)s'''

    P.run(statement)

###################################################
###################################################
###################################################

@merge(runSampleInference, "abundance.dir/merged_abundance.tsv")
def mergeAbundanceTables(infiles, outfile):
    '''
    combine sequence/abundance tables across
    samples
    '''
    statement = '''cgat combine_tables
                   -m 0
                   -c 1
                   -k 2
                   --use-file-prefix
                   -g abundance.dir/*_seq_abundance.tsv
                   --log=%(outfile)s.log
                   | sed 's/_seq_abundance.tsv//g'
                   > %(outfile)s'''
    P.run(statement)
    
#########################################
#########################################
#########################################

@follows(mkdir("taxonomy.dir"))
@transform(runSampleInference, regex("(\S+)/(\S+)_abundance.tsv"), r"taxonomy.dir/\2_taxonomy.tsv")
def assignTaxonomy(infile, outfile):
    '''
    assign taxonomy to sequences
    '''
    job_memory=PARAMS["taxonomy_memory"]
    taxonomy_file = PARAMS["taxonomy_taxonomy_file"]
    species_file = PARAMS["taxonomy_species_file"]

    statement = '''Rscript %(rscriptsdir)s/dada2_assign_taxonomy.R
                   --seqfile=%(infile)s
                   --training-set=%(taxonomy_file)s
                   --species-file=%(species_file)s
                   -o %(outfile)s'''
    P.run(statement)

###################################################
###################################################
###################################################

@transform(mergeAbundanceTables, suffix(".tsv"), "_id.tsv")
def addUniqueIdentifiers(infile, outfile):
    '''
    add unique identifiers to track the sequences
    '''
    outfile_map = outfile.replace(".tsv", ".map")
    PipelineDada2.seq2id(infile,
                         outfile_map,
                         outfile)
    
###################################################
###################################################
###################################################

@merge([assignTaxonomy, addUniqueIdentifiers],
       "taxonomy.dir/merged_taxonomy.tsv")
def mergeTaxonomyTables(infiles, outfile):
    '''
    combine sequence/taxonomy tables across
    samples
    '''
    PipelineDada2.mergeTaxonomyTables(infiles, outfile)
    

#########################################
#########################################
#########################################

@merge([mergeAbundanceTables, mergeTaxonomyTables],
       "abundance.dir/taxa_abundances.tsv")
def buildDefinitiveTable(infiles, outfile):
    '''
    build the final table with newids and
    abundance information
    '''
    PipelineDada2.makeDefinitiveAbundanceFile(infiles, outfile)

#########################################
#########################################
#########################################

@follows(mkdir("tree.dir"))
@transform(buildDefinitiveTable, regex("(\S+)/taxa_abundances.tsv"), r"tree.dir/tree.tsv")
def buildTree(infile, outfile):
    '''
    build a taxonomic tree of all ASVs that were detected across samples
    '''
    PipelineDada2.buildTree(infile, outfile)

#########################################
#########################################
#########################################

@follows(mkdir("taxonomy_abundances.dir"))
@split(buildDefinitiveTable, "taxonomy_abundances.dir/*_abundance.tsv")
def splitTableByTaxonomicLevels(infile, outfiles):
    '''
    split the table at different taxonomic levels and
    sum counts for each member of that level
    '''
    statement = '''Rscript %(rscriptsdir)s/dada2_split_levels.R
                   --outdir=taxonomy_abundances.dir
                   -i %(infile)s
                ''' 
    P.run(statement)

#########################################
#########################################
#########################################

@follows(buildTree, splitTableByTaxonomicLevels)
def full():
    pass


#########################################
#########################################
#########################################

# merge filter summaries into one table
@follows(filterAndTrim)
@merge(filterAndTrim, "filtered.dir/merged_filter_summary.tsv")
def mergeFilterSummary(infiles, outfile):
    '''
    Merging all filter summaries into one file
    reads.in    reads.out    sample
    '''
    
    PipelineDada2.mergeFilterSummary(infiles, outfile)


#########################################
#########################################
#########################################

# merge qc summaries into one table
@follows(runSampleInference)
@merge(runSampleInference, "abundance.dir/merged_qc_summary.tsv")
def mergeQCSummary(infiles, outfile):
    '''
    Merging all qc summaries into one file
    denoised    nochim    sample
    '''

    PipelineDada2.mergeQCSummary(infiles, outfile)

#########################################
#########################################
#########################################

# save pipeline yml files as tsv
@originate("parameter_table.tsv")
def yml2Table(outfile):
    
    ''' 
    Save parameters specified by pipeline.yml as tsv
    '''
    PipelineDada2.yml2Table(PARAMS, outfile) 

#########################################
#########################################
#########################################

# export tables into sqlite database for app compatibility
@transform([addUniqueIdentifiers, mergeTaxonomyTables,
            mergeFilterSummary, mergeQCSummary, yml2Table],
           regex(r"(.*)\.tsv"), r"\1.load")
def build_db(infile, outfile):
    '''
    Stores data generated throughout pipeline as a sqlite database.
    Structure of data tables and database is meant for compatibility
    with the shiny app
    '''
    
    # record merged_filter_summary, merged_qc_summary,
    # merged_taxonomy, merged_abundance_id
    # and yml table in database
    db_name = PARAMS["database_name"]
    P.load(infile, outfile, options=f'--database-url=sqlite:///./{db_name}')

#########################################
#########################################
#########################################

@follows(mkdir("report.dir"))
def build_report():
    '''
    render the rmarkdown report file
    '''
    reportdir = os.path.join(os.path.dirname(PipelineDada2.__file__), "pipeline_docs/Rmd/pipeline_dada2/")

    # cp a file as error model - random one really
    errF_files = glob.glob("filtered.dir/*errF.png")
    errF = errF_files[0]
    statement = '''cp %(errF)s report.dir/errF.png'''
    P.run(statement)
    
    # copy files to report directory
    statement = '''cd report.dir; cp %(reportdir)s/*.Rmd .; cd ../'''
    P.run(statement)

    # render the report
    statement = '''cd report.dir; R -e "rmarkdown::render('report.Rmd', output_file='report.html')"; cd ../'''
    P.run(statement)

    if PARAMS["report_diagnostics"] == 1:
        # render the diagnostic report
        statement = '''cd report.dir; R -e rmarkdown::render('pipeline_dada2_diagnostics.Rmd', output_file='pipeline_dada2_diagnostics.html'); cd ../'''
        P.run(statement)
    else:
        E.info("Not running diagnostics report")
        
#########################################
#########################################
#########################################
def main(argv=None):
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
