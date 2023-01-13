#!/usr/bin/env python
# coding: utf-8

from itertools import zip_longest

__version__ = "0.6"  

DNA_bases = "ATGC"
RNA_bases = "AUGC"

base_dict = {"A":"T","T":"A","G":"C","C":"G", "N":"N"}

def process(record):
    '''Creates a fastq record from a given list'''
    names = ["header", "seq", "+", "quality"]
    return {key:value for key, value in zip(names, record)}

def reverse_compl(seq):
    '''Returns reverse complement sequence of DNA strand'''
    rev_seq = ""    
    for base in seq: 
        rev_seq = base_dict[base] + rev_seq

    return (rev_seq)

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return(ord(letter)-33)

def qual_score(phred_score: str) -> float:
    '''This function calculates the average quality score of a given phred string'''
    assert(len(phred_score) != 0), "Phred schore is 0" #check to see if phred score is not 0

    sum_Qscore=0
    for letter in phred_score:
        sum_Qscore = sum_Qscore + convert_phred(letter)
    
    return (sum_Qscore/len(phred_score))

def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def validate_DNA_seq(DNA): 
    '''This function takes a string and returns True if string is a DNA sequence.'''
    DNA == DNA.upper() # makes everything uppercase 
    
    # if the base counts doesn't equal to the length, return false
    return len(DNA) == DNA.count("A") + DNA.count("T") + DNA.count("G") + DNA.count("C")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_DNA_seq(DNA), "String contains invalid characters"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def oneline_fasta(filename: str) -> str:
    '''# Takes in a fasta file and concatenate sequence lines into one 
    line for each fasta record. Outputs new file fasta and returns filename. '''
    
    fasta_dict = {} # dictionary to hold fasta header and single line sequence
    
    # parse through fasta file and strip new lines
    # concatenate sequence lines into one line
    with open(filename, "r") as file:
        sequence = ""
        header = ""

        for line in file:
            line = line.strip()

            if line.startswith(">"):
                header = line.strip()
                sequence =""

            else: 
                sequence += line.strip()
                fasta_dict[header] = sequence # update dictionary
                    
    # write out new one line sequence fasta file
    with open("oneLineFasta.fa", 'w') as file: 
        for header in fasta_dict:
            file.write(header + "\n" + fasta_dict[header] + "\n")

        return (file)
        
if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"

    assert validate_DNA_seq("AATAGAT") == True, "DNA string not recognized"
    assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"

    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
           
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
           
    assert qual_score("FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@") == 37.62105263157895, "wrong average phred score"

