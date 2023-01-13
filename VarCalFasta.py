#!/usr/bin/env python

# Version 0.1
# Compatible with Python Version 3.5

# This script reformats an input fasta file for UMI-VarCal.

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="input FASTA file")
    parser.add_argument("-o", "--output", required=True, type=str, help="output FASTA file")
    return(parser.parse_args())

args = get_args()
input_f = args.input
output_f = args.output

with open(input_f, "r") as fasta, \
open(output_f + ".fasta", "wt") as out:

    for line in fasta:
    
        line = line.strip()
      
        if line.strip().startswith(">"):
            line = line.split()[0][1:]
            print(">chr" + line, file=out)
        else:
           print(line, file=out)
    