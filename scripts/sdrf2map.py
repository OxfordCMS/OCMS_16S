'''
sdrf2map.py
====================================

:Author: Nick Ilott
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

An SDR file may contain metadata that relates to sample mapping. If the data come from 
the OGC then the file names contain the sample names. A KgenID column must also be present
for the script to run.

Usage
-----

.. Example use case

Example::

   python sdrf2map.py

Type::

   python sdrf2map.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import cgatcore.experiment as E
import glob

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--sdrf", dest="sdrf", type="string",
                      help="supply sdrf file")

    parser.add_option("-n", "--nlanes", dest="nlanes", type="int",
                      help="number of lanes sequenced")

    parser.add_option("-p", "--paired", dest="paired", action="store_true",
                      help="number of lanes sequenced")

    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    assert options.nlanes, "must specify the number of lanes to take ids from"
    
    sdrf = open(options.sdrf)
    header = sdrf.readline()[:-1].split("\t")
    assert "KgenID" in header, "sdrf file must contain kgenID column"

    # get the index of the kgenID
    kgenid_idx = [x for x, i in enumerate(header) if i == "KgenID"][0]
    
    for line in sdrf.readlines():
        data = line[:-1].split("\t")
        kgenid = data[kgenid_idx]

        # we use the number of lanes argument to take that many columns
        # from the end of the sdrf file - this is where they should be.
        # need to add a check for this. *2 for the 2 reads in a pair
        if options.paired:
            files = data[-(options.nlanes*2):]
            # don't need both reads in pairs for sample ids
            files = files[0::2]
        else:
            files = data[-(options.nlanes/2):]

        ids = ["_".join(x.split("_")[0:3]) for x in files]
        
        for i in ids:
            options.stdout.write("\t".join([i, kgenid]) + "\n")
        
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
