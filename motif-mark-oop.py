#!/usr/bin/env python

import argparse
import re

def parse_args():
    """
    Parses arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, help="name of input fasta containg sequences to be \
                        searched for motifs. Sequences must be <=1000 bp. Introns lower case, \
                        exons upper case. One exon and flanking introns")
    parser.add_argument('-m', type=str, help="name of input file containg motifs. One per line. \
                        ambiguous nucleotides permitted as denoted by IUPAC degenerate base system, \
                        upper or lower case, Ts and Us are both permissable\
                        allowed")
    return parser.parse_args()

def convert_to_regex(motif):
    """
    Takes in a motif string with ambiguous bases, returns string with regex meta characters to 
    account for IUPAC degenerate base system

    String will be all lower case, as the search function within the Sequence class searches by 
    converting all string within text to lowercase 
    """
    meta_motif = motif.lower()
    # convert any Us (or us) to Ts
    meta_motif = re.sub("u", "t", meta_motif)
    # w is a or t 
    meta_motif = re.sub("w", "[at]", meta_motif)
    # s is c or g
    meta_motif = re.sub("s", "[cg]", meta_motif)
    # m is a or c
    meta_motif = re.sub("m", "[ac]", meta_motif)
    # k is g or t
    meta_motif = re.sub("k", "[gt]", meta_motif)
    # r is a or g 
    meta_motif = re.sub("r", "[ag]", meta_motif)
    # y is c or t 
    meta_motif = re.sub("y", "[ct]", meta_motif)
    # b is c, g, or t 
    meta_motif = re.sub("b", "[cgt]", meta_motif)
    # d is a, g, t
    meta_motif = re.sub("d", "[agt]", meta_motif)
    # h is a, c, t
    meta_motif = re.sub("h", "[act]", meta_motif)
    # v is a, c, g
    meta_motif = re.sub("v", "[acg]", meta_motif)
    # n is all bases
    meta_motif = re.sub("n", "[actg]", meta_motif)
    return meta_motif


class Sequence():
    """
    Sequence to be searched, function to search
    """
    def __init__(self, sequence, exon_start_, exon_stop_, motifs_list):
        # The sequence itself
        self.seq =  sequence
        # 0-based indexing of the start of the exon - first capital letter
        self.exon_start = exon_start_
        # 0-based indexing of the end of the exon - first lowercase letter
        self.exon_stop = exon_stop_
        # Spot where each motif is found is stored as a dict.
        #   keys - the motif, as written in the txt file
        #   values - a list containing the start index (0-based) of where the motif was found in the sequnce
        motif_dict = dict()
        for motif in motifs_list:
            motif_dict[motif] = [] # intialize to empty list
        self.found = motif_dict

    def search(self):
        """
        Searches the sequence for each of the motifs which comprise motif_dict and updates the dict with the starting 
        index of each 
        """
        # Iterate over motifs
        for motif in self.found:
             # convert to all lower case 
            lower = self.seq.lower()
            # convert to regex expression
            re_motif = convert_to_regex(motif)
            # search lowercase string for each motif
            # start by getting number of motifs in string
            num_matches = len(re.findall(re_motif, lower))
            # string will be truncated continuouly, last keeps track of absolute index
            last = 0
            # iterating ocer matches for a given motif
            for i in range(num_matches):
                # Note that this strategy doesn't account for a single motif overlapping with itself
                # If that implementation is needed, a sliding window should be used, iterating over the string
                # one base at a time. This implementation also probably shouldn't use regex at all, given that 
                # regex matches can't overlap with one another but seeing that regex  was highly recommended 
                # for this assignment, I'm assuming it's cool.
                match = re.search(re_motif, lower)
                index_found = match.span()[0] + last
                self.found[motif].append(index_found)
                # add ending index to last to keep track of absolute index
                last += match.span()[1]
                # truncate string to get only chars after first match
                lower = lower[match.span()[1]:]


# Create object for the whole drawing, pass in all Input Sequences


def parse_fasta(filename):
    """
    Reads in input fasta, outputs sequences 
    """
    seq = ""
    seqs = []
    with open(filename) as fh:
        for line in fh:
            line=line.strip()
            if line[0]==">" and len(seq)>0:
                # everytime a new header line is reached, add seq to list and clear it
                seqs.append(seq)
                seq = ""
            elif not line[0]==">":
                seq = seq + line
    # add last sequence to list
    seqs.append(seq)
    return seqs

def read_motifs(filename):
    """
    Reads in each line, stores in list
    """
    mots = []
    with open(filename) as fh:
        for line in fh:
            line=line.strip()
            mots.append(line)
    return mots

# parse args
args = parse_args()

# Read in the files 
seqs = parse_fasta(args.f)
mots = read_motifs(args.m)

# initalize sequence objects for each sequence, store in list 
seq_objects = []
for seq in seqs:
    print(seq)
    # get starting and ending indices of exon
    match = re.search("[a-z][A-Z]",seq)
    # first capital letter is end of span
    start = match.span()[1]
    match = re.search("[A-Z][a-z]",seq)
    # first lowercase letter is end of span
    end = match.span()[1]
    seq_instance = Sequence(exon_start_=start, exon_stop_=end, motifs_list=mots, sequence=seq)
    seq_objects.append(seq_instance)

# search all of them for all motifs
for seq in seq_objects:
    seq.search()


