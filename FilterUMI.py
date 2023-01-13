#!/usr/bin/env python

# Read in 2 fastq files: R1 and R2

# Need to trim UMI/LINDA adapters from paired end reads.
    # Store these reads as headers in header line 
# Also need to find reverse complement UMI/LINDAs for reverse reads (R2)
    # Trim from R2 if found 

# Write out LINDA reads that we didn't find to a new file.
    # Need to create a read 1 and read 2 file
    # What to name these headers? Nothing? 

# Keep count in dictionary of how many UMI/LINDA combos found. 
    # Keep count of unidentified UMI/LINDAs


import numpy as np
import argparse
import re
import gzip
from itertools import zip_longest


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-R1", "--read1", required=True, type=str, help="read1 fastq filename")
    parser.add_argument("-R2", "--read2", required=True, type=str, help="read2 fastq filename")
    parser.add_argument("-e", "--end", required=True, type=int, help="first number of bps to search for umi")
    
    return(parser.parse_args())

args = get_args()
R1 = args.read1
R2 = args.read2
end = args.end

# create a set of LINDA sequences
linda_set = {"ACATGGCTGA", "GTTCAAGCAC", "TGGACCTACT", "CACGTTAGTG"}

# create dictionary for complement bases
#base_dict = {"A":"T","T":"A","G":"C","C":"G", "N":"N"}

# dictionary of count of each UMI/LINDA read found including unidentified LINDAs  
count_dict = {"ACATGGCTGA":0,"GTTCAAGCAC":0, "TGGACCTACT":0, "CACGTTAGTG":0,"unknown":0}

# set of unique UMI and LINDA combinations found 
umi_tag_set = set()

def process(record):
    names = ["header", "seq", "+", "quality"]
    return {key:value for key, value in zip(names, record)}

# Check if LINDA exists in R1 sequence line. If yes, trim UMI and LINDA from R1
# Return R1 record with trimmed sequence and UMI/LINDA in header if applicable
# and R2 record with updated header 
# Return True if LINDA sequences were found 
def find_UL(record1, record2):

    line = record1["seq"][0:end]
    for linda in linda_set:
        
        # found linda. find index start and end position
        # find UMI 8 bps before linda then trim starting at linda end position
        # trim quality score line at linda end position ,s
        # update count dictionary and return updated R1/R2 record

        if linda in line:
            linda_start = line.index(linda)
            linda_end = linda_start + 10

            if len(line[0:linda_start]) >= 8:
                umi = line[linda_start-8:linda_start]
                
                record1["seq"] = record1["seq"][linda_end:]
                record1["quality"] = record1["quality"][linda_end:]
                record1["header"] = f"{record1['header']}_{umi}_{linda}"
                record2["header"] = f"{record2['header']}_{umi}_{linda}"

                count_dict[linda] += 1
                if (umi+linda) not in umi_tag_set:
                    umi_tag_set.add(umi+linda)

                return record1, record2, True

    count_dict["unknown"] += 1
    return record1, record2, False

# return reverse complement of a given DNA string 
# def reverse_compl(seq):
#     rev_seq = ""    
#     for base in seq: 
#         rev_seq = base_dict[base] + rev_seq

#     return (rev_seq)


# write out fastq record to given file
def write_out(record, output):
    for key in record:
        print(record[key], file=output)

# write out all unique UMI/LINDA combinations found 
def write_umitag():
    with open("umi_tag_results.txt", "w") as out:
        
        print(f"{len(umi_tag_set)} UMIs found:", file=out)
        for umi in umi_tag_set:
            print(f"{umi}", file=out)

# write out a report of all LINDAs and unknown found
def calc_report():   
    with open("result.txt", "w") as out:
        num_records = sum(count_dict.values())
    
        # round(umi_count[key]/n * 100, 2)
        for key in count_dict: 
            percent = str(round(count_dict[key]/num_records * 100, 2))
            print(f"{key}:\t{count_dict[key]} reads\t{percent}%", file=out)            

# read R1 and R2 fastq files
# open output files

with gzip.open(R1, "rt") as R1_file,\
gzip.open(R2, "rt") as R2_file,\
gzip.open("UMI_trimmed_R1.fastq.gz", "wt") as out_R1,\
gzip.open("UMI_trimmed_R2.fastq.gz", "wt") as out_R2,\
gzip.open("unk_R1.fastq.gz", "wt") as unk_R1, \
gzip.open("unk_R2.fastq.gz", "wt") as unk_R2:

    # keep lists of current records
    record1 = []
    record2 = []
    
    for read1, read2 in zip_longest(R1_file, R2_file):

        record1.append(read1.strip())
        record2.append(read2.strip())

        # reached end of current records    
        if len(record2) == 4:
            record1 = process(record1)
            record2 = process(record2)
       
            # look for UMI and LINDA sequences in R1 
            record1, record2, found = find_UL(record1, record2)
        
            # if LINDAs weren't found, write out to unknown file
            if (found == True):
                write_out(record1, out_R1)
                write_out(record2, out_R2)
            else:
                write_out(record1, unk_R1)
                write_out(record2, unk_R2)

            # reset current record 
            record1 = []
            record2 = []    

calc_report()
write_umitag()