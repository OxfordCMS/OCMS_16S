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
import logging as L


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









#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
