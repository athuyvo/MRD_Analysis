#!/usr/bin/env python

import argparse
from matplotlib import pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", required=True, type=str, help="filepath to sam file")
    parser.add_argument("-p", "--probe", required=True, type=str, help="filepath to probe txt file")

    return(parser.parse_args())

args = get_args()
input_file = args.input
probe_file = args.probe

probe_dict = {}
no_probe = 0

# Create probe set with given probe text file 
def create_probe_set():
    with open(probe_file, "r") as f:
        for line in f:
            probe_dict[line.strip()] = 0

def write_count():
    with open("probe_counts.txt", "w") as f1:
        # for key in probe_dict:
        #     print(f"{key}:\t{probe_dict[key]}\tfound", file=f1)
        
        [print(f"{key}:\t{value}\tfound", file=f1) \
        for (key, value) in sorted(probe_dict.items(), key = lambda x:x[1])] 

        print(f"{sum(probe_dict.values())}\t found total", file=f1)
        print(f"{no_probe} no probe found", file=f1)      

with open(input_file, "r") as f, \
open("probe_filt.sam", "wt") as probe_f, \
open("no_probe.sam", "wt") as no_probe_f:
    create_probe_set()

    flag = ""
    for line in f:

        line = line.strip()    
        #print sam headers         
        if line.startswith("@") == True:
            print(f"{line}", file=probe_f)
            print(f"{line}", file=no_probe_f)
            continue

        seq = line.split("\t")[0:10]
        
        for probe in probe_dict:

            if probe in seq[9]:
                probe_dict[probe] +=1
                flag = "found"
        
        if flag == "found":
            print(f"{line}", file=probe_f)    
        else:
            print(f"{line}", file=no_probe_f)
            no_probe += 1

        flag = ""

write_count()

plt.bar(list(probe_dict.keys()), sorted(list(probe_dict.values())))
plt.title("Probe Distribution")
plt.xlabel("Probes")
plt.xticks([])
plt.ylabel("Num of Probes")
#plt.ylim(0,170000)
plt.savefig("probe.png")