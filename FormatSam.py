#!/usr/bin/env python

# Version 0.1

# This script reformats an input sam file and adds "chr" in front of chromosome number.

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="input sam file")
    parser.add_argument("-o", "--output", required=True, type=str, help="output sam file")
    return(parser.parse_args())

args = get_args()
input_f = args.input
output_f = args.output

with open(input_f, "r") as sam, \
open(output_f, "wt") as out:

    for line in sam:
            
        line = line.strip()

        # skip sam headers
        if line.startswith("@") == True:
            print(f"{line}", file=out)
            continue

        line = line.split("\t")
        line[2] = "chr" + line[2]
        print("\t".join(line), file=out)
            
        
        