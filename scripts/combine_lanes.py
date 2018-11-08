import sys
import os
import glob
import itertools

# input
directory = sys.argv[1]
wt2id = sys.argv[2]

def readWt2id(wt2id):
    '''
    read identifier map into dictionary
    '''
    mapping = {}
    inf = open(wt2id)
    for line in inf.readlines():
        data = line[:-1].split("\t")
        mapping[data[0]] = data[1]
    return(mapping)

def combineLanesRead1(directory, mapping_dict):
    '''
    read1 combination
    '''
    read1s = glob.glob(os.path.join(directory, "*_1.fastq.gz"))
    read1s.sort()
    found_index = set()
    out_map = open("read1.map", "w")
    for readfile in read1s:
        prefix = os.path.basename(readfile).replace("_1.fastq.gz", "")
        splitname = os.path.basename(readfile).split("_")
        lane, index = splitname[1], splitname[2]
        if index in found_index:
            continue
        else:
            found_index.add(index)
            tocombine = " ".join([x for x in read1s if index in x])
            newname = mapping_dict[prefix] + ".fastq.1.gz"
            out_map.write(tocombine + " " + newname + "\n")
            statement = f"""cat {tocombine} > {newname}"""
            os.system(statement)

def combineLanesRead2(directory, mapping_dict):
    '''
    read2 combination
    '''
    read2s = glob.glob(os.path.join(directory, "*_2.fastq.gz"))
    read2s.sort()
    found_index = set()
    out_map = open("read2.map", "w")
    for readfile in read2s:
        prefix = os.path.basename(readfile).replace("_2.fastq.gz", "")
        splitname = os.path.basename(readfile).split("_")
        lane, index = splitname[1], splitname[2]
        if index in found_index:
            continue
        else:
            found_index.add(index)
            tocombine = " ".join([x for x in read2s if index in x])
            newname = mapping_dict[prefix] + ".fastq.2.gz"
            out_map.write(tocombine + " " + newname + "\n")
            statement = f"""cat {tocombine} > {newname}"""
            os.system(statement)

wt2id = readWt2id(wt2id)
combineLanesRead1(directory, wt2id)
combineLanesRead2(directory, wt2id)
