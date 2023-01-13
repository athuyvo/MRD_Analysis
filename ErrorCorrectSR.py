#!/usr/bin/env python

# Version 0.1

# This is a genome referenced-based deduplication and error correction program. 
# Removes all PCR duplicates in a given sorted SAM file by Chromosome number then UMI. 
# This program for single-reads only. 



#################################################################
# Need to keep track of counts of each UMI. Remove less than 3. 
# Need to keep track of sequence to compare in each read family
# Remove any families that are less that 3 
#################################################################

##########################################################################
# UMI families need have be the same duplicates and be the same strand etc
##########################################################################

import argparse
import re



# set of tuples to store unique read info for current unique UMI: ({UMI, chrom, pos, strand})
# if record is in umipos_set, create family dictionary 
umipos_set = set()


# dictionary to keep track of records for current chromosome  
# {"UMI":None, "Chrom":None, "Pos":None, "Strand":None, "Seq":None, "Count": 0}
record_dict = {}


# unique reads per chromosome
count_dict = {"header":0, "uniq":0, "incorrect":0, "dups":0, "removed": 0}

# keep track of which chromosome we're checking.
curr_chrom = ""
prev_chrom = ""

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--input", required=True, type=str, help="absolute filepath to sorted .sam file")
    parser.add_argument("-o", "--output", required=True, type=str, help="absolute filepath to output sam file")
    parser.add_argument("-u", "--umi", required=True, type=str, help="absolute filepath to umi txt file")

    return(parser.parse_args())

args = get_args()
input_file = args.input
out_file = args.output
umi_file = args.umi


# Create UMI set with given UMI text file 
def create_UMI_set():
    with open(umi_file, "r") as f:
        line_count = 0
        for line in f:
            line_count += 1
            if line_count > 2 : 
                umi_set.add(line.strip())
    # umi_set.add("ACT")    # for test files 
    # umi_set.add("AAA")    # for test files

# parse bitwise flag
# check if sequence is reverse 
def check_flag(flag):

    # check if read is unmapped, or supplementary alignment
    # if yes, skip 
    # if ((flag & 4) == 4 or (flag & 2064) == 2064 or (flag & 2048) == 2048) :
    #     count_dict["removed"] += 1
    #     return "SKIP"


    # check if sequence is reverse
    # if yes, add "-" to current search list
    
    if ((flag & 16) == 16):
        record_dict["Strand"] = "-"
    else:
        record_dict["Strand"] = "+"
    
    # return ""

# parse sequence record for QNAME, bitwise flag, chromosome, position,
# and cigar string
# return "SKIP" to stop searching if UMI doesn't exist
def parse_sam(line):

    # Grab UMI in QNAME of sequence line (col 0)
    # Check if UMI is in UMI set
    # if yes, add to current search list
    # otherwise, ignore line and return "SKIP"

    # umi = re.search("[ATGC]{3}", line[0])[0] # for test files 
    umi = line[0].split(":")[-1]

    if umi in umi_set:
        record_dict["UMI"] = umi
    else: 
        count_dict["incorrect"] +=1
        return "SKIP"

    
    # Grab bitwise flag (col 1)
    # flag = check_flag(int(line[1]))
    check_flag(int(line[1]))
    # if flag == "SKIP":
    #     return "SKIP"    

    # chromosome RNAME (col 2)
    # add chromosome num to current search list
    record_dict["Chrom"]= line[2]
    prev_chrom = curr_chrom

    # add left-most position (col 3) to current search list 
    record_dict["Pos"] = int(line[3])

    # return cigar string (col 5)

    return line[5]


# parse input cigar string
# store alignment match, soft clipping, skipped regions, and deletion variables
# update position if needed 

def parse_cigar(cigar):

    # if current strand is positive, check for soft clipping on 5'
    # adjust left-most position if needed 

    if record_dict["Strand"] == "+":
        left_soft = [int(i) for i in re.findall("[0-9]+(?=S[0-9]+)", cigar)]
        record_dict["Pos"] -= sum(left_soft)


    # if strand is reverse, adjust for:
    # alignment match(es), deletions, right soft clippings, skipped regions
    # position = position + alignment match + skipped region + right alignment match 
    # + right soft clipping +  deletion - 1 (left adjusted, inclusive) ** didn't adjust for inclusive. 
    if record_dict["Strand"]== "-":
        record_dict["Pos"] += sum([int(i) for i in re.findall("[0-9]+(?=M)", cigar)]) # align
        record_dict["Pos"] += sum([int(i) for i in re.findall("[0-9]+(?=D)", cigar)]) # deletion
        record_dict["Pos"] += sum([int(i) for i in re.findall("(?<=[A-Z]{1})[0-9]+(?=S)", cigar)]) # right soft clip
        record_dict["Pos"] += sum([int(i) for i in re.findall("[0-9]+(?=N)", cigar)]) # skipped

# write out counter for headers, unique reads, duplicates, and incorrect UMIs found 
# write out counter for each unique read found per chromosome 
def write_counters():
    with open("deduper_report.txt", "w") as f1:
        for key in count_dict:
            print(f"{key}:\t{count_dict[key]}\tfound", file=f1)

# Read sorted input .sam file (will be sorted by chromosome then by UMI).
# Open write file. 

with open(input_file, "r") as f, \
open(out_file, "w") as out:
    create_UMI_set()
    for line in f:
    
        # If current line is a header, write out "@" headers to .sam file
        # else save line to write out when unique record is found 
        # update count dictionary with number of header lines
        
        if line.strip().startswith("@") == True:
            print(f"{line.strip()}", file=out)
            count_dict["header"] += 1 
            continue

        record = line.strip().split("\t")[0:9]
        cigar = parse_sam(record)

        # if UMI in sam record doesn't exist, stop search and continue to next record 
        # otherwise, continue parse through cigar string 
        if cigar == "SKIP":
            continue 
        parse_cigar(cigar)

        # add UMI, chromosome, updated position, and strand type to dictionary as key 
        # add sequence record as value and true for unique found
        # write out record to sam file if unique
        # update count dictionary of unique or duplicates 

        if tuple(record_dict.values()) not in umipos_set: 
           umipos_set.add(tuple(record_dict.values()))
           print(f"{line.strip()}", file=out)
           count_dict["uniq"] +=1         

        else: 
            # duplicate found
            # update duplicate counter 
            count_dict["dups"] +=1 

        # done searching for duplicates for this chromosome
        if prev_chrom != curr_chrom:
            #umipos_set = {}
            umipos_set.clear()
        
        record.clear()
    # done searching for all duplicates 
    # write out counters
    write_counters()