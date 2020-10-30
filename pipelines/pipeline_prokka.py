##############################################
##############################################
##############################################
# Simple pipeline to run prokka on multiple
# assemblies
##############################################
##############################################
##############################################

from ruffus import *
import cgatcore.pipeline as P
import cgatcore.iotools as IOTools
import os
import sys
import collections
import cgatcore.experiment as E

PARAMS = P.get_parameters(filenames=["pipeline.yml"])

##############################################
##############################################
##############################################

@follows(mkdir("annotations.dir"))
@transform("*.fna.gz", regex("(\S+).fna.gz"), r"annotations.dir/\1/\1.tsv")
def runProkka(infile, outfile):
    '''
    run prokka annotations - 
    very basic with no parameterisation
    '''
    newdirname  = os.path.join("annotations.dir", P.snip(infile, ".fna.gz"))
    job_memory = PARAMS["prokka_memory"]
    job_threads = PARAMS["prokka_threads"]
    prefix = P.snip(infile, ".fna.gz")
    t = P.snip(infile, ".gz")
    
    statement = '''zcat %(infile)s > %(t)s;
                   /gfs/devel/nilott/prokka/bin/prokka --cpus %(job_threads)s --outdir %(newdirname)s --prefix %(prefix)s %(t)s --force;
                   rm -rf *.fna
                '''
    P.run(statement)

##############################################
##############################################
##############################################

@merge(runProkka, "annotations.dir/merged_annotations.tsv")
def mergeAnnotations(infiles, outfile):
    '''
    merge the annotations into a matrix with a union of
    genes and presence/absence calls for each species
    '''

    # This might need to be optimised as was a
    # quick and dirty approach to aggregating data
    # Actually could use merge_tables from cgat?
    
    genes = set()
    for infile in infiles:
        inf = open(infile)
        # skip header
        inf.readline()
        for line in inf.readlines():
            data = line[:-1].split("\t")
            gene = data[3]
            if gene == "":
                continue

            # prokka annotates multiple copies with
            # _[0-9] so jsut take gene name for
            # comparison purposes. Note that this means
            # in the ouput that there will be no indication
            # as to whether multiple copies are present in
            # a genome.
            if len(gene.split("_")) == 2:
                gene = gene.split("_")[0]
            genes.add(gene)
        inf.close()

    # iterate over files again and add
    # presence/absence calls
    result = collections.defaultdict(list)
    colnames = list()
    outf = open(outfile, "w")

    for g in genes:
        for infile in infiles:
            # per file genes
            genomegenes = set()
            inf = open(infile)

            # skip header
            inf.readline()
            for line in inf.readlines():
                data = line[:-1].split("\t")
                gene = data[3]
                if len(gene.split("_")) == 2:
                    gene = gene.split("_")[0]
                genomegenes.add(gene)
            if g in genomegenes:
                result[g].append("1")
            else:
                result[g].append("0")

    outf.write("\t".join(["gene"] + [P.snip(os.path.basename(x), ".tsv") for x in infiles]) + "\n")
    for gene, data in result.items():
        outf.write("\t".join([gene] + data) + "\n")
    outf.close()
    
##############################################
##############################################
##############################################

@follows(mergeAnnotations)
def full():
    pass

##############################################
##############################################
##############################################

if __name__== "__main__":
    P.main()

    
