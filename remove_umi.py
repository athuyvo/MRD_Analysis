#!/usr/bin/env python
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="This program will take two paired-end fastq files and will remove UMIs and Linda sequences and add them to the header for deduping purposes")
    parser.add_argument("-r1", "--fastqr1", help="Input the filename for Read 1 (fastq)", type=str, required=True)
    parser.add_argument("-r2", "--fastqr2", help="Input the filename for Read 2 (fastq)", type=str, required=True)
    parser.add_argument("-or1", "--outputr1", help="Input the output filename for Read 1 (fastq)", type=str, required=True)
    parser.add_argument("-or2", "--outputr2", help="Input the output filename for Read 2 (fastq)", type=str, required=True)
    parser.add_argument("-of", "--outputfolder", help="Input an output folder for unknown reads, report and umi file", type=str, required=True)
    return parser.parse_args()

args = get_args()
fast = args.fastqr1
fast2 = args.fastqr2
out_r1_file = args.outputr1
out_r2_file = args.outputr2
o_folder = args.outputfolder

#VARIABLES
umi_count = {"ACATGGCTGA":0, "GTTCAAGCAC":0,"TGGACCTACT":0, "CACGTTAGTG":0} #Linda Sequence 25-28 respectively 
full_umi = set()
not_umi = 0
no_linda_file_r1 = o_folder + "unknown_r1.fastq.gz"
no_linda_file_r2 = o_folder + "unknown_r2.fastq.gz"
report= o_folder + "umi_report.txt"
umi_file= o_folder + "known_umi.txt"


#FUNCTIONS

def read_file(f:str) -> list:
    '''This function will take in a fastq file and return a list containing the data of a record'''
    head = f.readline().strip()
    seq = f.readline().strip()
    sep = f.readline().strip()
    qsc = f.readline().strip()
    if head == "":
        return 0
    return [head,seq,sep,qsc]


def find_umi_r1(read:str, line:int) -> list:
    '''This function will take in a record and revise if any of the linda sequences is on the selected range and if any it will return an umi(with linda) and the sequence without the umi'''
    umi = ""
    new_seq = ""
    #for every linda in linda dictionary
    for linda in umi_count.keys():
        #if there is a linda in the range of 0:21
        if linda in read[8:25]:
            #Will only cut the first linda found
            x = read.split(linda,1)
            umi = x[0] + linda
            new_seq = x[1]
            if len(umi) > 18:
                umi = umi[-18:]
                new_seq = x[1]
                seq_len=len(new_seq)
                umi_count[linda] += 1
            elif len(umi) == 18:
                umi = umi
                new_seq = x[1]
                seq_len=len(new_seq)
                umi_count[linda] += 1
            elif len(umi) < 18:
                return(False, None, None, None)
            if umi not in full_umi:
                full_umi.add(umi)
            return(True, umi, new_seq, seq_len)
            break
    else:
        return(False, None, None, None)


def header_change(head:str, umi_c:str, qscr:str, length:int) -> str:
    '''It will read the header and the corresponding linda sequence and will return a new header and correct the quality score'''
    #new_head = head + ":" + umi_c #split head by white space, append to [0], join back the list with an underscore before the umi
    split = head.split()
    #new_head = split[0] + ":" + umi_c + " " + split[1]
    new_head = split[0] + "_" + umi_c + " " + split[1]
    if length != len(qscr):
        new_qscr = qscr[:length]
    else: 
        new_qscr=qscr
    return(new_head, new_qscr)

def write_out(linda_exist:bool, no_lindaf:str, normalf:str, record:list, umi:str, seq:str, length:int):
    '''This function will take in a True or False statement and depending on which, will write out to a no_linda file or a new fastq file with header, sequence and quality score updated'''
    if linda_exist == False:
            no_lindaf.write("\n".join(record)+"\n")
    else:
        record[0], record[3] = header_change(record[0], umi, record[3], length)
        record[1] = seq
        normalf.write("\n".join(record)+"\n")

#Open all the files that are going to be written into
f_nlr1 = gzip.open(no_linda_file_r1,"wt")
f_nlr2 = gzip.open(no_linda_file_r2,"wt")
f_or1 = gzip.open(out_r1_file, "wt")
f_or2 = gzip.open(out_r2_file, "wt")


#ACTUAL CODE
#Open two fastq files
n=0
with gzip.open(fast,"rt") as fhr1, gzip.open(fast2,"rt") as fhr2:
    while True:
        R1=read_file(fhr1)
        R2=read_file(fhr2)
        #If the read 1 encounters an empty line it will stop reading and break the loop
        if R1 == 0:
            break
        n+=1
        #head,sequence,separator,quality_scr = R1
        #Obtain if the record R1 has a linda, an umi and the new sequence
        linda_umi_seq_len = find_umi_r1(R1[1], n)
        if linda_umi_seq_len[0] == False:
            not_umi += 1
        
        #write out with the arguments. Theres a linda?, file if no linda, file if linda, Original record, UMI, New Sequence
        write_out(linda_umi_seq_len[0], f_nlr1, f_or1, R1, linda_umi_seq_len[1], linda_umi_seq_len[2], linda_umi_seq_len[3])
        #Now for R2
        #write out with the arguments. Theres a linda?, file if no linda, file if linda, Original record, UMI of R1, original Sequence
        write_out(linda_umi_seq_len[0], f_nlr2, f_or2, R2, linda_umi_seq_len[1], R2[1], len(R2[1]))
    

with open(report,"w") as fr:
    fr.write("Stats from Removing UMIs\n\n")
    fr.write("The Ammount of reads is: " + str(n) + "\n\n")
    for key in umi_count:
        fr.write(str(key) + "\t" + str(umi_count[key]) + "\t" + str(round(umi_count[key]/n * 100, 2)) + "%\n\n")
    fr.write("The ammount of records without full UMIs are:" + str(not_umi))

with open(umi_file,"w") as fu:
    fu.write("List of UMIs detected\n\n")
    for thing in full_umi:
        fu.write(thing + "\n")
            

f_nlr1.close()
f_or1.close()
f_nlr2.close()
f_or2.close()




#THIS THAT DIDN'T MAKE IT IN THE CODE BUT MIGHT BE VALUABLE FOR THE FUTURE

#THIS FUNCTION WON'T BE NECESSARY AS WE WON'T BE LOOKING FOR THE UMITAG IN R2
# def find_umi_r2(read2:str) -> list:
#     umi = ""
#     for linda in umi_count_r2.keys():
#         if linda in read2[-20:]:
#             x = read2.split(linda)
#             if umi == "": #For the first line wont have a value
#                 #Define which will the UMI be
#                 umi = linda + x[1]
#                 umi = rev_complement(umi)
#                 linda = rev_complement(linda)
#                 new_seq = x[0]
#                 if umi not in full_umi:
#                     full_umi.add(umi)
#                     umi_count[linda] += 1
#             return(True, umi, new_seq)
#             break
#     else:
#         return(False, None, None)

#DON'T NEED TO REVERSE COMPLEMENT ANYTHING

# rev = {"A": "T", "C": "G", "T": "A", "G":"C", "N" : "N"}
# def rev_complement(seq: str) -> str:
#     """This function will read a sequence and return it's reverse complement"""
#     rev_com = ""
#     for letter in seq:
#         rev_com += rev[letter]
#     return(rev_com[::-1])

#THIS IS TO WRITE OUT TO THE FILE, BUT I MADE A FUNCTION THAT DOES THIS

# if linda_umi_seq[0] == False:
#     f_nl.write("\n".join(R1)+"\n")
# else:
#     R1[0], R1[3] = header_change(R1[0], linda_umi_seq[1], R1[3])
#     R1[1] = linda_umi_seq[2]
#     f_or1.write("\n".join(R1)+"\n")