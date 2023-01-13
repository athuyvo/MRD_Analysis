#!/usr/bin/env python

# Reverse complements a given text file. 

import bioinfo
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="input file")
    parser.add_argument("-o", "--output", required=True, type=str, help="output file")
    return(parser.parse_args())

args = get_args()
input_f = args.input
output_f = args.output

with open(input_f, "r") as f, \
open(output_f, "wt") as out:

    for line in f:
        seq = line.strip()
        newline = bioinfo.reverse_compl(seq)
        print(f"{newline}", file=out)
