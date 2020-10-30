################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
pipeline_subsamplebam
===========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline to randomly subsample a set of bamfiles using the lowest count out of all
samples.

Uses samtools to achieve this.

Code
====

"""
from ruffus import *

import sys
import glob
import collections
import re
import os
import time
import itertools
import optparse, shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import numpy as np

P.getParameters( 
    [ "pipeline.ini" ] )

PARAMS = P.PARAMS

################################################
################################################
################################################

@follows(mkdir("readcounts.dir"))
@transform("*.bam",
           regex("(\S+).bam"),
           r"readcounts.dir/\1.readcount")
def countMappedReads(infile, outfile):
    '''
    count the number of mapped reads per sample
    '''
    statement = '''samtools view -F4 %(infile)s | wc -l > %(outfile)s'''
    P.run()

################################################
################################################
################################################
    
@follows(mkdir("subsampled.dir"))
@transform("*.bam",
           regex("(\S+).bam"),
           add_inputs(countMappedReads),
           r"subsampled.dir/\1.subsampled.bam")    
def subsampleBam(infiles, outfile):
    '''
    subsample bam file based on lowest counts. This is
    converted to a fraction for subsampling using
    samtools
    '''
    bamfile = infiles[0]
    
    # match counts file with bam file
    track = P.snip(infiles[0], ".bam")
    counts_files = infiles[1:]
    counts_file = [c for c in counts_files if track in c][0]

    # get the minimum number of mapped reads
    mapped_counts = []
    for inf in counts_files:
        count = open(inf).readline().strip()
        mapped_counts.append(count)
    min_count = np.amin([x for x in map(int, mapped_counts)])

    # get fraction to subsample
    count = open(counts_file).readline().strip()
    fraction = np.round(float(min_count)/float(count),2)

    # subsample with samtools
    statement = '''samtools view -s %(fraction)f -b %(bamfile)s > %(outfile)s;
                   checkpoint;
                   samtools index %(outfile)s'''
    P.run()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )    
