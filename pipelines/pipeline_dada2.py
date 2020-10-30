'''
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
import pipelines.PipelineDada2 as PipelineDada2

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# Initialize the pipeline
P.initialize()
PARAMS = P.get_params()

# scripts directory - R scripts for dada2 functions
scriptsdir = os.path.basename(pipeline_dada2.__file__).replace("pipelines", "R")

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
        statement = '''zcat %(inf)s > %(outtmp)s'''
        P.run(statement)

    # hackabout
    outtmp = os.path.join(tmpdir, [x for x in infiles if x.endswith(".fastq.1.gz")][0].replace(".gz", ""))
    statement = '''Rscript %(scriptsdir)s/dada2_filter_and_trim.R
                           --infile=%(outtmp)s
                           %(paired)s
                           --maxN=%(maxn)s
                           --maxEE=%(maxee)s
                           --truncQ=%(truncQ)s
                           --truncLen=%(truncLen)s
                           --trimLeft=%(trimLeft)s
                           --filtered-directory=filtered.dir'''
    P.run(statement)
    summary = open(outfile.replace(".fastq.1.gz" + "_summary.tsv")).readline()
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

    if PARAMS["paired"] == 1:
        # get second read in pair
        infile_read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")

        # generate statement to execute
        statement = '''Rscript %(scriptsdir)s/dada2_sample_inference.R
                       --filtF=%(infile)s
                       --filtR=%(infile_read2)s
                       --nbases=%(nbases)s
                       %(options)s
                       --outdir=%(outdir)s'''
    else:
        statement = '''Rscript %(scriptsdir)s/dada2_sample_inference.R
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

    statement = '''Rscript %(scriptsdir)s/dada2_assign_taxonomy.R
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
    statement = '''Rscript %(scriptsdir)s/dada2_split_levels.R
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
def build_db(infiles, outfile):
    '''
    Stores data generated throughout pipeline as a sqlite database.
    Structure of data tables and database is meant for compatibility
    with the shiny app
    '''
    
    # record merged_filter_summary, merged_qc_summary,
    # merged_taxonomy, merged_abundance_id
    # and yml table in database
    P.load(infiles, outfile)

    # add index to each table
    # statement = '''sqlite3 csvdb; CREATE UNIQUE INDEX merged_taxonomy ON merged_abundance (sequence); .quit'''
    # P.run(statement)

#########################################
#########################################
#########################################

@follows(mkdir("report.dir"))
def build_report():
    '''
    render the rmarkdown report file
    '''
    reportdir = os.path.basename(pipeline_dada2.__file__).replace("pipelines", "docs/Rmd/pipeline_dada2")
    author = '"' + PARAMS["report_author"] + '"'
    title = '"' + PARAMS["report_title"] + '"'

    # cp a file as error model - random one really
    errF_files = glob.glob("filtered.dir/*errF.png")
    errF = errF_files[0]
    statement = '''cp %(errF)s report.dir/errF.png'''
    P.run(statement)
    
    # copy files to report directory
    statement = '''cd report.dir; cp %(reportdir)s/*.Rmd .; cd ../'''
    P.run(statement)

    # add report_author and report_title
    statement = '''sed -i 's/report_author/%(author)s/g' report.dir/report.Rmd;
                   sed -i 's/report_title/%(title)s/g' report.dir/report.Rmd'''
    P.run(statement)

    # add report_title
    statement = '''sed -i 's/report_title/%(title)s/g' report.dir/pipeline_dada2.Rmd'''
    P.run(statement)

    # render the report
    statement = '''cd report.dir; R -e "rmarkdown::render(report.Rmd, output_file=report.html)"; cd ../'''
    P.run(statement)

    if PARAMS["report_diagnostics"] == 1:
        # render the diagnostic report
        statement = '''cd report.dir; R -e rmarkdown::render(pipeline_dada2_diagnostics.Rmd, output_file=pipeline_dada2_diagnostics.html); cd ../'''
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
