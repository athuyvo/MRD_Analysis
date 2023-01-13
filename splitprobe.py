#!/usr/bin/env python

# This script is used to remove the common 5' sequence from MRD probes. 

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="input probe file")
    parser.add_argument("-o", "--output", required=True, type=str, help="output probe file")
    return(parser.parse_args())

args = get_args()
input_f = args.input
output_f = args.output

with open(input_f, "r") as fi, \
open(output_f, "wt") as out:

    for line in fi:
    
        line = line.strip()
        line = line.split("ATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
        line = line[1]
        print(line, file=out)
    