#!/usr/bin/env python

# Author: Anh Vo 
# University of Oregon 
# Version 2.0

# This script is designed for Illumina Paired-End and Single reads. It takes an aligned sam file and 
# removes reads that are not mapped in proper pairs and removes supplementary alignnments. 
# Reads that have not been filtered out are written to a new sam file.

# It takes into account reads that may not be in proper paired order in paired-end sam files and 
# the possibility that there may proper reads that are be a missing its pair. In this case, 
# these reads will be written out to an unmatched sam file. 

# Reads that do not pass the above criteria will be written out to an unmatched sam file. 

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--input", required=True, type=str, help="input sam file")
    parser.add_argument("-pe", "--pe", required=True, choices=["True", "False"], help="if read is PE or SR")

    return(parser.parse_args())


args = get_args()
input_file = args.input
pe = args.pe

# if pe is not True and pe is not False:
#     exit("FileTypeError: Must indicate -pe as True or False")


# Returns True PE reads are properly paired
def check_proper(flag):
    if ((flag & 2) == 2):
        return True
    return False

# Returns True if reads are not supplementary alignments
def check_supp(flag):
    if ((flag & 2048) != 2048):
        return True 
    return False

# Returns True if reads are aligned
def check_unmapped(flag):
    if ((flag & 4) != 4): 
        return True
    return False 

# Filter through PE reads in sam file for mapped, properly paired, and 
# removes supplementary alignments
# Reads that matched the criterias are written out to a filtered
# sam file else written to unmatched sam file
def run_pe():

    # dict to keep track of records to write out 
    r1_proper= {}
    r2_proper= {}

    with open(input_file , "r") as f, \
    open("FilteredAligned.sam", "w") as out:
        for line in f: 
            
            record = line.strip()

            # skip sam headers
            if line.startswith("@") == True:
                print(f"{record}", file=out)
                continue
            
            flag = int(line.strip().split("\t")[1])
            cluster = split_record(record)
            
            # Find R1 
            # Save to list to write out if proper paired and not supp read
            # Otherwise, write out to unmatched sam file
            if ((flag & 64) == 64):  # R1
                #proper = check_proper(flag)

                # if ((flag & 2) == 2 and (flag & 2048) != 2048):
                if check_proper(flag) and check_supp(flag):
                
                    # Check R2 dictionary in case R2 proper paired is out of order in sam file
                    # Else add to R1 dictionary
                    if cluster in r2_proper:
                        print(f"{record}\n{r2_proper[cluster]}", file=out)
                        del r2_proper[cluster]
                    else:
                        r1_proper[cluster] = record
                continue

            # Found R2
            # If R2 mate (R1) is properly paired, print out to proper matched file
            if ((flag & 128) == 128): # R2
                if cluster in r1_proper:
                    print(f"{r1_proper[cluster]}\n{record}", file=out)
                    del r1_proper[cluster]

                # Check in case R1 proper paired is out of order in sam file 
                else:
                    if check_proper(flag) and check_supp(flag):
                    # if ((flag & 2) == 2 and (flag & 2048) != 2048):
                        r2_proper[cluster] = record

    write_unmatched(r1_proper)
    write_unmatched(r2_proper)

# Filter through single read sam file and remove 
# unmapped and supplementary reads
def run_sr():
    with open(input_file , "r") as f, \
    open("FilteredAligned.sam", "w") as out:
    #open("unmatched.sam", "w") as unm_out:

        for line in f: 
            record = line.strip()

            # skip sam headers
            if line.startswith("@") == True:
                print(f"{record}", file=out)
                continue
            
            flag = int(line.strip().split("\t")[1])

            if check_supp(flag) and check_unmapped(flag): 
                print(f"{record}", file=out)

# Split qname into clusters to help identify R1 and R2 pairs 
def split_record(line):
    qname = line.strip().split("\t")[0]
    cluster_x = qname.split(":")[5] 
    cluster_y = qname.split(":")[6]
    cluster_y = cluster_y.split("_")[0]
    
    return cluster_x + ":" + cluster_y
    

# Write out remaining that reads that were not 
# matched to its proper pair to an unmatched sam file
def write_unmatched(dct):
    with open("unmatched.sam", "w") as unm_out:
        for value in dct.values():
            print(f"{value}", file=unm_out)


if pe == "True":
    run_pe()

if pe == "False":
    run_sr()
                    
