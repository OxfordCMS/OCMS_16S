"""
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

The input is a directory of fastq formatted files.

"""

# load modules
from ruffus import *

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import shutil
import logging as L
import os
import sys
import glob

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["pipeline.ini"])

PARAMS = P.PARAMS

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

    # make sure parameters are present
    assert maxn != "" and maxee != "" and truncQ != "" and truncLen != "", \
    "must specify all parameters to filterAndTrim"
    
    if PARAMS["paired"] == 1:
        paired = "--paired"
    else:
        paired = ""

    tmpdir = P.getTempDir()

    # unzip one by one
    infile_read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")
    infiles = [infile, infile_read2]
    
    for inf in infiles:
        outtmp = inf.replace(".gz", "")
        statement = '''zcat %(inf)s > %(tmpdir)s/%(outtmp)s'''
        P.run()
        
    statement = '''Rscript %(scriptsdir)s/dada2_filter_and_trim.R
                   --infile=%(infile)s
                   %(paired)s
                   --maxN=%(maxn)s
                   --maxEE=%(maxee)s
                   --truncQ=%(truncQ)s
                   --truncLen=%(truncLen)s
                   --directory=%(tmpdir)s
                   --filtered-directory=filtered.dir'''
    P.run()
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
    nreads = PARAMS["sample_inference_nreads"]
    outdir = os.path.dirname(outfile)

    if PARAMS["paired"] == 1:
        # get second read in pair
        infile_read2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")

        # generate statement to execute
        statement = '''Rscript %(scriptsdir)s/dada2_sample_inference.R
                       --filtF=%(infile)s
                       --filtR=%(infile_read2)s
                       --nreads=%(nreads)s
                       --outdir=%(outdir)s'''
    else:
        statement = '''Rscript %(scriptsdir)s/dada2_sample_inference.R
                       --filtF=%(infile)s
                       --nreads=%(nreads)s
                       --outdir=%(outdir)s'''

    P.run()

###################################################
###################################################
###################################################

@merge(runSampleInference, "abundance.dir/merged_abundance.tsv")
def mergeAbundanceTables(infiles, outfile):
    '''
    combine seqquence/abundance tables across
    samples
    '''
    statement = '''python %(cgatscriptsdir)s/combine_tables.py
                   -m 0
                   -c 1
                   -k 2
                   --use-file-prefix
                   -g abundance.dir/*_seq_abundance.tsv
                   --log=%(outfile)s.log
                   | sed 's/seq_abundance.tsv//g'
                   > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
