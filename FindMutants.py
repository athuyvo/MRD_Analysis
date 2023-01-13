#!/usr/bin/env python

# This script will find known tumor mutants in a SAM or FASTQ using a provide sequence txt file. 
# This will work for aligned reads or SE/MERGED FASTQ files only.
# Output a fastq or sam file with double mutants founds

import argparse
from itertools import zip_longest
import statistics as stat

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-m", "--mutant", required=True, type=str, help="mutant sequence file")
    parser.add_argument("-f", "--file", required=True, type=str, help="input sam file")
    parser.add_argument("-t", "--type", required=True, type=str, help="fastq or sam file")

    return(parser.parse_args())

args = get_args()
mutant_file = args.mutant
input_file = args.file
file_type = args.type

# FILE MUST BE IN SAM OR FASTQ FORMAT OR EXIT PROGRAM
input_filetype = input_file.split(".")[1]                  

if file_type != "fastq" and file_type != "sam" and file_type != "fq":
    exit("FileTypeError: File must be in fastq or sam format")

if input_filetype != file_type:
    exit("FileTypeError: Input filetype doesn't not match specified filetype")

# chrom: [wt, dbl, pos49, pos 52]
# dict containing double mutants at position 49 and 42 and mutant chromosome position 
wt_dict = {}
mut_dict = {} 
pos49_dict = {}
pos52_dict = {}

# dict containing summary of mutant counts: {mutant type: Mean count, Mean %,  Median count, Median %}
sum_dict = {"WT": [None, None, None], "Double": [None, None, None], "Pos 49": [None, None, None], "Pos 52": [None, None, None]} 

# tuple count of single mutants and each mutant chromosome position 
wt_count = {}
mut_count = {}
pos49_count = {}
pos52_count = {}



# Create mutant sets for single mutants 
# Create double mutant dictionary with chromosome position as value and sequence as key
def create_mutant_set():

    with open(mutant_file, "r") as f:
        skip = True

        for line in f:    

            # skip header line
            if skip == True:
                skip = False
                continue 

            line = line.strip().split("\t")

            # wt_dict[line[0]] = line[1]
            mut_dict[line[0]] = [line[1], line[2], line[3], line[4]]
            #pos49_dict[line[0]] = line[3]
            #pos52_dict[line[0]] = line[4]
            wt_count[line[0]] = 0
            mut_count[line[0]] = 0
            pos49_count[line[0]] = 0
            pos52_count[line[0]] = 0 

        



# find mutants in sequence 
def find_mutant(seq):
    flag = False
#    print(f"Seq: {seq}")
    global mut_dict
    global wt_count
    global mut_count
    global pos49_count
    global pos52_count
    
    for mut in mut_dict:
        for index in range(0,4):

        # search for double mutants
        ##### fix this ###### 
            if mut_dict[mut][index] in seq: 
                if index == 0: 
                    wt_count[mut] += 1
                
                # found double mutant
                if index == 1:
                    mut_count[mut] += 1
                    flag = True

                if index == 2:
                    pos49_count[mut] += 1 
                
                if index == 3: 
                    pos52_count[mut] += 1 
        
    return flag

def calc_stats(mut, dct):  
    sum_dict[mut][0] = sum(dct)
    sum_dict[mut][1] = round(stat.mean(dct),4)
    sum_dict[mut][2] = stat.median(dct)



# write out total mutant counts 
# write out only chromosome position counts for double mutants
def write_count(num_records):
    with open("mutants_found.tsv", "wt") as f1:
        # chrom: [wt, dbl, pos49, pos 52]
        # sum_dict = sum, mean, median
        # sum_dict = {"WT": None, None, None, "Double": None, None, None, "Pos 49": None, None, None, "Pos 52": None, None, None } 


        # sum up all values in mutant and wild type counts 
       
        calc_stats("WT", wt_count.values())
        calc_stats("Double", mut_count.values())
        calc_stats("Pos 49", pos49_count.values())
        calc_stats("Pos 52", pos52_count.values())
        

    
        ###### fix the sum and count here
        
        print(f"Chrom Position\tNum Wild Type\t% Wild Type\tNum Double Mutants\t% Double Mutants\tNum Pos 49\t% sPos 49\tNum Pos 52\t% Pos 52", file = f1)

        # for key in zip(wt_count, mut_count, pos49_count, pos52_count):
        for key in mut_dict:
            # percent of mutants found is divided by number of recordsfound for that chromosome position            

            wt_percent = 0 
            mut_percent = 0
            pos49_percent = 0
            pos52_percent = 0

            total = wt_count[key] + mut_count[key] + pos49_count[key] + pos52_count[key]

            if total > 0 :
                wt_percent = round(wt_count[key] / total * 100, 4)
                mut_percent = round(mut_count[key] / total * 100, 4)
                pos49_percent = round(pos49_count[key] / total * 100, 4)
                pos52_percent = round(pos52_count[key] / total * 100, 4)
            
    
            print(f"{key}:\t{wt_count[key]}\t{wt_percent}%\t{mut_count[key]}\t{mut_percent}%\t{pos49_count[key]}\t{pos49_percent}%\t{pos52_count[key]}\t{pos52_percent}%", file = f1)

        print(f"Total Mutants:\t{sum_dict['Double'][0]}\tMean:\t{sum_dict['Double'][1]}\tMedian:\t{sum_dict['Double'][2]}", file = f1)
        print(f"Total Pos 49:\t{sum_dict['Pos 49'][0]}\tMean:\t{sum_dict['Pos 49'][1]}\tMedian:\t{sum_dict['Pos 49'][2]}", file = f1)
        print(f"Total Pos 52:\t{sum_dict['Pos 52'][0]}\tMean:\t{sum_dict['Pos 52'][1]}\tMedian:\t{sum_dict['Pos 52'][2]}", file = f1)
        print(f"Total WT:\t{sum_dict['WT'][0]}\tMean:\t{sum_dict['WT'][1]}\tMedian:\t{sum_dict['WT'][2]}", file = f1)
        print(f"Total number of records:\t{num_records}", file=f1)


# find mutants in sam file 
def run_sam(f):
    with open("mutants.sam", "wt") as out:
        num_records = 0
        
        for line in f:
            
            line = line.strip()

            # skip sam head
            if line.startswith("@") == True:
                print(f"{line}", file=out)
                continue

            num_records += 1
            seq = line.split("\t")[9]
            
            flag = find_mutant(seq)

            # if double mutant found, print out record
            if flag == True:
                print(line, file=out)

    return num_records


def process(record):
    names = ["header", "seq", "+", "quality"]
    return {key:value for key, value in zip(names, record)}

# find mutant in fastq file
def run_fastq(f):

    with open("mutants.fastq", "wt") as out:
        num_records = 0
        curr_record = []

        for line in f:

            curr_record.append(line.strip())
   
            # end of current record 
            if len(curr_record) == 4:
                num_records += 1
                record = process(curr_record)
                curr_record = []

                # search for mutant
                flag = find_mutant(record["seq"])
                
                # if double mutant found, print out record
                if flag == True:
                    for key in record:
                        print(record[key], file=out)
            
    return num_records

with open(input_file, "r") as f:
    create_mutant_set()
    search = ""
    flag = False
    
    if file_type == "sam":
       num_records = run_sam(f)
       write_count(num_records)
    
    if file_type == "fastq" or file_type == "fq":
        num_records = run_fastq(f)
        write_count(num_records)



