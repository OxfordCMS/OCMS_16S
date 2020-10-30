'''
combine_featurecounts_summaries.py 
====================================

:Author:
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Simple script to combine featurecounts summaries.

Usage
-----

.. Example use case

Example::

   python combine_featurecounts_summaries.py

Type::

   python combine_featurecounts_summaries.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import cgatcore.experiment as E
import glob


def getTotalCount(infile):
    '''
    get total count from file for outputting
    proportions
    '''
    inf = open(infile)
    inf.readline()
    total = 0
    for line in inf.readlines():
        data = line[:-1].split("\t")
        total += int(data[1])
    inf.close()
    return(total)


def reformat(infile):
    '''
    Not sure yet
    '''
    total = getTotalCount(infile)
    inf = open(infile)
    sample = os.path.basename(inf.readline()[:-1].split("\t")[1])
    for line in inf.readlines():
        data = line[:-1].split("\t")
        proportion = (float(data[1])/total)*100
        sys.stdout.write("\t".join(data + [str(proportion), sample]) + "\n")
    inf.close()
        
def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--directory", dest="directory", type="string",
                      help="supply directory where the input summaries aer located")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    infiles = glob.glob(os.path.join(options.directory, "*/*genes*summary*"))
    sys.stdout.write("category\tnreads\tpreads\tsample\n")
    for infile in infiles:
        reformat(infile)

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
