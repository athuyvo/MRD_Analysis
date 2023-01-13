#!/usr/bin/env python

# Version 0.1
# Compatible with Python Version 3.5

# This script formats an input BED file for UMI-VarCal. 

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="input BAM file")
    parser.add_argument("-o", "--output", required=True, type=str, help="output BAM file")
    return(parser.parse_args())

args = get_args()
input_f = args.input
output_f = args.output

with open(input_f, "r") as bed, \
open(output_f + ".bed", "wt") as out:

    for line in bed:
        line = line.strip()
        chrom = "chr" + line.split()[0]
        print(chrom, line[1:], sep="", file=out)
