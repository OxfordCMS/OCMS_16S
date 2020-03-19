"""
=======================================================
Filter fastq files based on reads unmapped in bam file
=======================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

"""

# load modules
from ruffus import *

import cgatcore.experiment as E
import logging as L

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random
import cgat.FastaIterator as FastaIterator
import cgat.Fastq as Fastq

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import cgatcore.pipeline as P
P.get_parameters(
    ["pipeline.ini"])

PARAMS = P.PARAMS

###################################################################
###################################################################
###################################################################


@follows(mkdir("filtered.dir"))
@transform(glob.glob("*.fastq.gz"), 
           regex("(\S+).fastq.gz"),
           add_inputs(glob.glob("*.bam")),
           r"filtered.dir/filtered.\1.fastq.gz")
def filterFastq(infiles, outfile):
    '''
    filter fastq file based on read tyep 
    speicifed in bam file
    '''
    t = PARAMS.get("reads_type")
    fastq = infiles[0]
    track = P.snip(fastq, ".fastq.gz")
    bams = infiles[1]
    if len(bams) == 1:
        bam = bams[0]
    else:
        bam = [bam for bam in bams if bam.find(track) != -1][0]

    job_options="-l mem_free=30G"
    
    if PARAMS.get("reads_invert") == 1:
        invert = "--invert"
    else:
        invert = ""
    statement = '''zcat %(fastq)s | python %(scriptsdir)s/fastq2filteredfastq.py
                                    -b %(bam)s
                                    %(invert)s
                                    --reads=%(t)s
                                    --log=%(outfile)s.log
                                    | gzip > %(outfile)s'''
    P.run()

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
