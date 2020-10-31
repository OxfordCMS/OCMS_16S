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
import pandas as pd

def seq2id(seqtable, outfile_map, outfile_table):
    '''
    assign a unique identifier to each unique
    sequence in dada2 output
    '''
    inf = open(seqtable)
    header = inf.readline().strip('\n').split('\t')
    header = ["featureID"] + header[1:]
    c = 1
    out_map = open(outfile_map, "w")
    out_map.write("id\tfeatureID\n")
    out_table = open(outfile_table, "w")
    out_table.write('\t'.join(header) + "\n")
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

def makeDefinitiveAbundanceFile(infiles, outfile):
    '''
    make a definitive table of abundances of the form:
    ASVXXX:k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species    21
    ...
    '''
    # read in taxonomy file (contains all id keys)
    tax_file = open(infiles[1], 'r')
    tax_file.readline()

    # read in sequence counts (id by sequence)
    seq_file = open(infiles[0], 'r')
    seq_file.readline()

    # create dictionary for identifiers
    id_dict = {}

    for line in tax_file.readlines():
        data = line.strip("\n").split("\t")
        values = [x for i, x in enumerate(data) if i != 1]
        id_dict[data[1]] = values
    tax_file.close()


    # initiate output file with modified header from merged_abundance
    outfile = open(outfile, "w")
    with open(infiles[0], 'r') as f:
        header = f.readline().strip("\n").split("\t")
        header[0] = 'featureID'
    outfile.write("\t".join(header) + "\n")
    
    # iterate over abundance file and create new names
    for line in seq_file.readlines():
        data = line.strip("\n").split("\t")

        # look up sequence in id dictionary made from
        # merged taxonomy file
        identifier = data[0]
        id_value = id_dict.get(identifier)
        newid = id_value[0] + ":" + id_value[-1]
    
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
    
    pattern = r".*(?=_seq_abundance.tsv)"
    result = [m.group(0) for x in fname for m in [re.search(pattern, x)] if m]
    summary_file = [x + "_summary.tsv" for x in result]

    # get header from first summary file
    first_file = open(os.path.join("abundance.dir",
                                   summary_file[0]),'r')
    header = first_file.readline()
    # initiate merged summary file
    merged = open(outfile, "w")
    merged.write(header)

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
# save pipeline.yml file as dataframe (tsv)

def yml2Table(param_dict, outfile):

    # expand parameters dictionary so only one value per key
    # keep record parameter sections
    key_list = []
    val_list = []
    task_list = []
    
    for key in param_dict.keys():
        
        curr_val = param_dict[key]

        # if value is a dictionary, unpack dictionary values
        if isinstance(curr_val, dict):
            
            for k in curr_val.keys():
                v = curr_val[k]

                # convert value to list of strings
                ## when nested dict value is one value
                ## immediately update output lists
                if isinstance(v, list) == False:
                    key_list.append(k)
                    val_list.append(str(v))
                    task_list.append(key)
                ## when nested dict value is list
                ## split up list so one value per key
                else:
                    for i, x in enumerate(v):
                        entry_key = str(k) + str(i)
                        key_list.append(entry_key)
                        val_list.append(x)
                        task_list.append(key)

        # assuming only single values in non-dictionary values
        # so can add directly to output lists
        else:
            key_list.append(str(key))
            val_list.append(str(curr_val))
            task_list.append(key)
    
    # record pipeline.yml into a table
    yml_table = pd.DataFrame(list(zip(task_list, key_list, val_list)),
                             columns = ["task", "parameter","value"])

    # write table to file
    yml_table.to_csv(outfile, index = False, sep = "\t")
   
