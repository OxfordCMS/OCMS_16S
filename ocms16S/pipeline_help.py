'''
pipeline_help.py
====================================================

:Author: Nick Ilott
:Tags: Python

Purpose
-------

run -h for the dada2 rscripts. Provides options for use in pipeline.yml.

Usage
-----

.. Example use case

Example::

   ocms_16s help --filterAndTrim

Type::

   ocms_16s help --filterAndTrim

for command line help.

Command line options
--------------------

'''

import ocms16S.PipelineDada2 as PipelineDada2
import sys
import os
import cgatcore.experiment as E


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--filterAndTrim", dest="filter_and_trim",
                        action='store_true', help="return options for filterAndTrim task")
    parser.add_argument("--sampleInference", dest="sample_inference",
                        action='store_true', help="return options for sampleInference task")
    parser.add_argument("--assignTaxonomy", dest="assign_taxonomy",
                        action='store_true', help="return options for assignTaxonomy task")

    # add common options (-h/--help, ...) and parse command line
    (args) = parser.parse_args()

    rscriptsdir = os.path.join(os.path.dirname(PipelineDada2.__file__), "R")
    if args.filter_and_trim:
        assert not args.sample_inference and not args.assign_taxonomy, "can only view options for one task at a time"
        script = os.path.join(rscriptsdir, "dada2_filter_and_trim.R")
    elif args.sample_inference:
        assert not args.filter_and_trim and not args.assign_taxonomy, "can only view options for one task at a time"
        script = os.path.join(rscriptsdir, "dada2_sample_inference.R")
    elif args.assign_taxonomy:
        assert not args.sample_inference and not args.filter_and_trim, "can only view options for one task at a time"
        script = os.path.join(rscriptsdir, "dada2_sample_inference.R")
    os.system(f'Rscript {script} -h')


if __name__ == "__main__":
    sys.exit(main(sys.argv))
