'''
cms_basic.py
=============

:Author: Nick Ilott
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Takes directories of where basic OCMS data analyses were performed and creates a new directory that contains output files that are to
be transferred to the customer.

Files that are provided are:

Fastqc report (html)
dada2 report (html)
taxa_abundances.tsv (ASV level count data)
species_abundances.tsv (count data at species level)
genus_abundances.tsv (count data at genus level)
family_abundances.tsv (count data at family level)
order_abundances.tsv (count data at order level)
class_abundances.tsv (count data at class level)
phylum_abundances.tsv (count data at phylum level)


Usage
-----

.. Example use case

Example::

   python cms_basic.py --dada2-dir=<path-to-dada2-run> --fastqc-dir=<path-to-fastqc-run> --project-name=<name-of-project>

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

def buildFileDescription(project_name):
    '''
    build a file description file
    '''
    file_layout=f"""
    |-{project_name}
    |------Reports
    |---------dada2_report.html
    |---------fastqc_report.html
    |------Data
    |---------taxa_abundances.tsv
    |---------species_abundances.tsv
    |---------genus_abundances.tsv
    |---------family_abundances.tsv
    |---------order_abundances.tsv
    |---------class_abundances.tsv
    |---------phylum_abundances.tsv
    |---------merged_abundance_id.map

    dada2_report.html -> Metrics on the dada2 run including reads input and output and taxonomic assignments.
    fastqc_report.html -> Quality summary of raw sequencing data.
    taxa_abundances.tsv -> Read counts for each amplicon sequence variant (ASV)(rows) detected per sample (columns).
    species_abundances.tsv -> Read counts for each species detected i.e. ASVs collapsed at the level of species annotation.
    genus_abundances.tsv -> Read counts for each genus detected i.e. ASVs collapsed at the level of genus annotation.
    family_abundances.tsv -> Read counts for each family detected i.e. ASVs collapsed at the level of family annotation.
    order_abundances.tsv -> Read counts for each order detected i.e. ASVs collapsed at the level of order annotation.
    class_abundances.tsv -> Read counts for each class detected i.e. ASVs collapsed at the level of class annotation.
    phylum_abundances.tsv -> Read counts for each class detected i.e. ASVs collapsed at the level of phylum annotation.
    merged_abundance_id.map -> Text file that maps the ASV identifier (arbitrary number) to the original identified sequence byt dada2.
    """
    outf = open(f"{project_name}/Files.txt", "w")
    outf.write(file_layout)
    outf.close()

############################################
############################################
############################################
    
def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--dada2-dir", dest="dada2_dir", type="string",
                      help="supply dada2 run directory")
    parser.add_option("--fastqc-dir", dest="fastqc_dir", type="string",
                      help="supply fastqc run directory")
    parser.add_option("--project-name", dest="project_name", type="string",
                      help="project name used to create master directory")

    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    dada2_dir = options.dada2_dir
    fastqc_dir = options.fastqc_dir
    project_name = options.project_name
    
    # make directories
    os.system(f"mkdir {project_name}; mkdir {project_name}/Reports; mkdir {project_name}/data")

    # get all of the dada2 data
    dada2_data = [os.path.join(dada2_dir, "abundance.dir/taxa_abundances.tsv"), os.path.join(dada2_dir, "abundance.dir/merged_abundance_id.map")] 
    dada2_data = dada2_data + glob.glob(os.path.join(dada2_dir, "taxonomy_abundances.dir/*"))
    dada2_data = " ".join(dada2_data)
    dada2_report = os.path.join(dada2_dir, "report.dir/report.html")

    # do the copying
    os.system(f"cp {dada2_data} {project_name}/data; cp {dada2_report} {project_name}/Reports; mv {project_name}/Reports/report.html {project_name}/Reports/dada2_report.html")

    # copy the fastq file over
    os.system(f"cp {fastqc_dir}/MultiQC_report.dir/multiqc_report.html {project_name}/Reports; mv {project_name}/Reports/multiqc_report.html {project_name}/Reports/fastqc_report.html")

    # build file description
    buildFileDescription(project_name)

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
