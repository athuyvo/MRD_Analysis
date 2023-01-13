#!/usr/bin/env python

# Author: Anh Vo 
# University of Oregon 
# Version 1.0

# This script is designed for Illumina Paired-End reads. It takes an aligned sam file and 
# removes reads that are not mapped in proper pairs and removes supplementary alignnments. 
# Reads that have not been filtered out are written to a new sam file.

# It takes into account reads that may not be in proper paired order in the sam file and 
# the possibility that there may proper reads that are be a missing its pair. 

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--input", required=True, type=str, help="input sam file")

    return(parser.parse_args())


args = get_args()
input_file = args.input

# dict to keep track of records to write out 
r1_proper= {}
r2_proper= {}

def check_proper(flag):
    if ((flag & 2) == 2 and (flag & 2048) != 2048):
        return True
    return False

# Write out remaining that reads that were not 
# matched to its proper pair to an unmatched sam file
def write_unmatched(dct, out_file):
        for value in dct.values():
            print(f"{value}", file=out_file)

# Read input sam file and write out properly paired reads to new sam file     
with open(input_file , "r") as f, \
open("FilteredAligned.sam", "w") as out:

    # Check for R1 or R2 for each line 
    # If R1, filter out nonproper paired reads and supplementary reads 
    # If R2, remove if R1 was removed
    for line in f: 
        record = line.strip()
        
        # Split qname into clusters to help identify R1 and R2 pairs 
        qname = line.strip().split("\t")[0]
        cluster_x = qname.split(":")[5] 
        cluster_y = qname.split(":")[6]
        cluster_y = cluster_y.split("_")[0]
        cluster = cluster_x + ":" + cluster_y 


        flag = int(line.strip().split("\t")[1])
        
        # Found R1 
        # Save to list to write out if proper paired and not supp read
        # Otherwise, write out to unmatched sam file
        if ((flag & 64) == 64):  # R1
            proper = check_proper(flag)
            #if ((flag & 2) == 2 and (flag & 2048) != 2048):
            if proper:
            
                # Check R2 dictionary in case R2 proper paired is out of order in sam file
                # Else add to R1 dictionary
                if cluster in r2_proper:
                    print(f"{record}\n{r2_proper[cluster]}", file=out)
                    del r2_proper[cluster]
                else:
                    r1_proper[cluster] = record
            continue

            # else: 
            #     r1_unm[cluster] = record
            #     current = cluster
            #     # find r2 to remove and write out to unmatched file
            #     continue

        # Found R2
        # If R2 mate (R1) is properly paired, print out to proper matched file
        if ((flag & 128) == 128): # R2
            if cluster in r1_proper:
                print(f"{r1_proper[cluster]}\n{record}", file=out)
                del r1_proper[cluster]

            # Check in case R1 proper paired is out of order in sam file 
            else:
                if ((flag & 2) == 2 and (flag & 2048) != 2048):
                    r2_proper[cluster] = record

with open("unmatched.sam", "w") as unm_out:
    write_unmatched(r1_proper, unm_out)
    write_unmatched(r2_proper, unm_out)
    print(r1_proper)
    print(r2_proper)
                    


#check if r1 == 64
 #
# remove:


# supplementary reads == 2048
# not proper paired != 2

# if removing: 
# check next read for r2 == 128 and remove too
#check qname for same UMI tag or x and y cluster. should be the same for both 
     # if not, save into dictionary? and rrevisit or find mate?
     # for now, print out to unmatched file 





# remove singletons and umi family of 2 or less (later)
