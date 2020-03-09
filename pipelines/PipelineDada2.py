##################################################
##################################################
##################################################
# functions for use with pipeline_dada2.py
##################################################
##################################################
##################################################

import re
import os
import itertools

def seq2id(seqtable, outfile_map, outfile_table):
    '''
    assign a unique identifier to each unique
    sequence in dada2 output
    '''
    inf = open(seqtable)
    header = inf.readline()
    c = 1
    out_map = open(outfile_map, "w")
    out_map.write("id\tsequence\n")
    out_table = open(outfile_table, "w")
    out_table.write(header)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        seq = data[0]
        identifier = "ASV" + str(c)
        c += 1
        out_map.write("\t".join([identifier, seq]) + "\n")
        out_table.write("\t".join([identifier, "\t".join(data[1:])]) + "\n")
    out_map.close()
    out_table.close()
        
##################################################
##################################################
##################################################

def mergeTaxonomyTables(infiles, outfile):
    '''
    merge taxonomy files
    '''

    # read in mapping file
    id_dict = {}
    id_map = open("abundance.dir/merged_abundance_id.map", 'r')
    id_map.readline()
    for line in id_map.readlines():
        entry = line.strip("\n").split("\t")
        id_dict[entry[1]] = entry[0]

    # taxonomy files
    pattern = r"seq_taxonomy\.tsv"
    tax_files = [x for x in infiles if re.search(pattern, x)]

    taxonomy = {}
    for infile in tax_files:
        inf = open(infile, "r")
        inf.readline()
        for line in inf.readlines():
            data = line[:-1].split("\t")
            seq, tax = data[0], [data[1], data[2], data[3], data[4], data[5], data[6], data[7]]

            # strip quotes from any taxon names
            tax = [x.strip('"') for x in tax]

            # replace missing data with NA
            tax = ["NA" if x=='' else x for x in tax]
            
            # keep joined version of taxa to add as Taxon column; add on tax level prefix
            prefix = ['k','p','c','o','f','g','s']
            taxon = [prefix[i] + "__" + x for i,x in enumerate(tax)]
            taxon = ";".join(taxon)

            tax.append(taxon)

            # remove the __ in names - will allow consistency
            # downstream
            # tax = [re.sub(r"\S__", "", x) for x in tax]

            # check there are the correct number of elements
            assert len(tax) == 8, "Not enough levels in %s" % ", ".join(tax)
            taxonomy[seq] = tax

        inf.close()

    # write merged taxonomy file
    outf = open(outfile, "w")
    outf.write("featureID\tsequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tTaxon\n")
    for seq, taxa in taxonomy.items():
        entry = id_dict[seq] + "\t" + seq + "\t" + "\t".join(taxa) + "\n"
        outf.write(entry)
    outf.close()


##################################################
##################################################
##################################################

def makeDefinitiveAbundanceFile(id2seq, seq2taxonomy, id2abundance, outfile):
    '''
    make a definitive table of abundances of the form:
    ASVXXX:p__phylum;c__class;o__order;f__family;g__genus;s__species    21
    ...
    '''
    inf_id2seq = open(id2seq)
    inf_id2seq.readline()
    inf_seq2taxonomy = open(seq2taxonomy)
    inf_seq2taxonomy.readline()
    inf_id2abundance = open(id2abundance)
    header = inf_id2abundance.readline()
    
    id2seq_dict = {}
    seq2taxonomy_dict = {}

    for line in inf_id2seq.readlines():
        data = line[:-1].split("\t")
        id2seq_dict[data[0]] = data[1]

    for line in inf_seq2taxonomy.readlines():
        data = line[:-1].split("\t")
        # from phylum level
        taxonomy = data[2:]
        taxonomy[0] = "p__" + taxonomy[0]
        taxonomy[1] = "c__" + taxonomy[1]
        taxonomy[2] = "o__" + taxonomy[2]
        taxonomy[3] = "f__" + taxonomy[3]
        taxonomy[4] = "g__" + taxonomy[4]
        taxonomy[5] = "s__" + taxonomy[5]
        taxonomy = ";".join(taxonomy)
        
        seq2taxonomy_dict[data[0]] = taxonomy

    # iterate over abundance file and create
    # new names
    outfile = open(outfile, "w")
    outfile.write(header)
    for line in inf_id2abundance.readlines():
        data = line[:-1].split("\t")
        identifier = data[0]
        seq = id2seq_dict[identifier]
        taxonomy = seq2taxonomy_dict[seq]
        newid = ":".join([identifier, taxonomy])
        outfile.write("\t".join([newid] + data[1:]) + "\n")
    outfile.close()
    

##################################################
##################################################
##################################################

def buildTree(infile, outfile):
    '''
    builds a text format tree file for putting into
    graphlan. Takes merged_taxonomy.tsv as input
    just go down to genus - Removes 
    '''
    inf = open(infile)
    inf.readline()

    outf = open(outfile, "w")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        asv_taxon = data[0]
        asv_taxon = asv_taxon.replace(":", ";")
        asv_taxon = asv_taxon.split(";")
        asv_taxon = asv_taxon[1:6] + [asv_taxon[0]]
        asv_taxon = ".".join(asv_taxon)
        outf.write(f'{asv_taxon}\n')
    outf.close()

##################################################
##################################################
##################################################

def mergeFilterSummary(infiles, outfile):

    '''
    merge all filter summaries in filter.dir
    into merge_filter_summary.tsv
    '''

    # identify filter summaries
    fname = [os.path.basename(x) for x in infiles]
    pattern = r".*(?=\.fastq)"
    result = [m.group(0) for x in fname for m in [re.search(pattern, x)] if m]
    summary_file = [x + "_summary.tsv" for x in result]

    # initiate merged summary file
    merged = open(outfile, "w")
    merged.write("reads.in\treads.out\tsample\n")

    # read one summary file at a time
    for f in summary_file:
        with open(os.path.join('filtered.dir', f), "r") as curr_file:
            entry = curr_file.readlines()[1]
            merged.write(entry)
        curr_file.close()
    merged.close()

###################################################
###################################################
###################################################
def mergeQCSummary(infiles, outfile):
    '''
    merge all qc summaries into one file
    into merge_qc_summary.tsv)
    '''
    
    # identify filter summaries
    fname = [os.path.basename(x) for x in infiles]
    print(fname)
    pattern = r".*(?=_seq_abundance.tsv)"
    result = [m.group(0) for x in fname for m in [re.search(pattern, x)] if m]
    summary_file = [x + "_summary.tsv" for x in result]
    # initiate merged summary file
    merged = open(outfile, "w")
    merged.write("denoised\tnochim\tsample\n")

    # read one summary file at a time
    for f in summary_file:
        with open(os.path.join('abundance.dir', f), "r") as curr_file:
            entry = curr_file.readlines()[1]
            merged.write(entry)
        curr_file.close()
    merged.close()

#####################################################
#####################################################
#####################################################
                     
