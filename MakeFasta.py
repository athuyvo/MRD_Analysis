#!/usr/bin/env python

# Version 0.1
# Compatible with Python Version 3.5


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

    count = 0
    for line in f:
        count += 1
        print(">",count, sep="", file=out)
        print(line.strip(), file=out)
    