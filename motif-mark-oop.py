#!/usr/bin/env python

import argparse
import cairo
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


def parse_fasta(filename):
    """
    Reads in input fasta, outputs sequences, names in list form
    """
    seq = ""
    seqs, names = [], []
    with open(filename) as fh:
        for line in fh:
            line=line.strip()
            if line[0]==">":
                # everytime a new header line is reached, add seq to list and clear it
                if len(seq)>0:
                    # Don't append empty seq on first iteration
                    seqs.append(seq)
                seq = ""
                names.append(line[1:])
            else:
                seq = seq + line
    # add last sequence to list
    seqs.append(seq)
    return seqs, names

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


class Sequence():
    """
    Sequence to be searched, function to search
    """
    def __init__(self, sequence, exon_start_, exon_stop_, motifs_list, name_):
        # The sequence itself
        self.seq =  sequence
        self.seq_len = len(sequence)
        self.name = name_
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
        # dict of motifs and overlaps, populated after search() is run. See refactor_found
        self.overlap = []

    def search(self):
        """
        Searches the sequence for each of the motifs which comprise motif_dict and updates the dict with the starting 
        index of each motif
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
    
    def refactor_found(self):
        """
        Once found dict is filled out, refactors the dict into a list of the following form:
        [start, end, motif, number motifs overlapped]
            start - start position (inclusive)
            end - end position (inclusive)
            motif - as written in txt file
            number motifs overlapped - 0 if no overlap, 1 if overlaps one other, 2 if overlaps 2 other, etc.
        """
        unsorted_list = []
        keys = self.found.keys()       
        for key in keys:
            if len(self.found[key])>0:
                # only make tuples out of motifs that have been found
                for start_pos in self.found[key]:
                    # populate with zero to start with
                    mini_list = [start_pos, start_pos+len(key)-1, key, 0]
                    unsorted_list.append(mini_list)

        # Now, sort all lists by start position
        sorted_list_of_lists = sorted(unsorted_list, key=lambda x: x[0])
        # two motifs are overlapping if the end of motif i >= to the start of motif i + 1
        # end_index is motif i, start_index is motif i+1
        end_index = 0
        start_index = 1
        while start_index<len(sorted_list_of_lists) and end_index<len(sorted_list_of_lists):
            # if the end of the first motif overlaps the start of the next
            if sorted_list_of_lists[end_index][1]>=sorted_list_of_lists[start_index][0]:
                # then increment the number of overlaps both have
                sorted_list_of_lists[end_index][3] += 1
                sorted_list_of_lists[start_index][3] += 1
                # move onto the next motif, continue incrementing this variable util the overlap
                # with motif i is cleared
                start_index += 1
            else:
                # if there is no overlap, move onto the next pair of motifs
                end_index += 1
                start_index = end_index + 1
        
        # store results in class
        self.overlap = sorted_list_of_lists

# Create object for the whole drawing, pass in all Input Sequences, figure name
class Drawing():

    def __init__(self, seqs, name_, motifs):
        self.seqs = seqs
        self.name = name_
        len_seq = 0
        # get longest sequence for drawing
        for seq in seqs:
            if seq.seq_len > len_seq:
                len_seq = seq.seq_len
        self.longest = len_seq
        self.motif_colors = dict()
        for index, motif in enumerate(motifs):
            # Max 5 motifs
            # starts at bluish magenta, goes to deep orange
            color = (203/255, 20/255, (255-index*50)/255)
            self.motif_colors[motif] = color

    def draw_seqs(self):
        """
        Draws introns and exons of sequence in black, label
        Draw motifs
        """
        edge = 50
        num_seq = len(self.seqs)
        width = self.longest + edge*2
        height_per_seq = 200
        height = num_seq*height_per_seq # excluding height needed for key
        exon_height = 50
        intron_height = 20

        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)
        # Create white background
        ctx.set_source_rgb(1, 1, 1)  
        ctx.paint()

        # draw exons and introns, label
        ctx.set_source_rgb(0,0,0) 
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(18)
        for i in range(num_seq):
            # draw height is upper left hand corner of intron
            draw_height = i*height_per_seq + height_per_seq/2 - intron_height
            # draw height - 50 is text position
            ctx.move_to(edge, draw_height-edge)
            ctx.show_text(self.seqs[i].name)
            intron_length = self.seqs[i].seq_len
            exon_length = self.seqs[i].exon_stop - self.seqs[i].exon_start
            # draw intron
            # (x, y, width, height) -> origin is TOP LEFT
            ctx.rectangle(edge, draw_height, intron_length, intron_height)
            ctx.fill() 
            ctx.rectangle(self.seqs[i].exon_start + edge, draw_height-(exon_height-intron_height)/2, exon_length, exon_height)
            ctx.fill()
        
        # Draw motifs
        for index, seq in enumerate(self.seqs):
            # i_draw_height is y coordinate of intron
            i_draw_height = index*height_per_seq + height_per_seq/2 - intron_height
            # e_draw_height is y coordinate of exon
            e_draw_height = i_draw_height-(exon_height-intron_height)/2
            # iterating over motif sequences by key - motif text
            prev_overlaps = 0 # number of consecutive overlaps happening before a given motif
            for motif in seq.overlap:
                # lists are (start, end, motif, # overlaps)
                start, end, motif_text, num_overlaps = motif
                # to show overlap, make motifs shorter by dividing by scaling factor
                scaling_factor = num_overlaps+1
                c1, c2, c3 = self.motif_colors[motif_text]
                ctx.set_source_rgb(c1, c2, c3) 
                # if motif comes before or after exon
                if (start < seq.exon_start and end < seq.exon_start) or (start > seq.exon_stop and end > seq.exon_stop):
                    # draw motif intron width
                    motif_height = intron_height/scaling_factor
                    m_draw_height = i_draw_height
                else:
                    # otherwise it's drawn exon width
                    motif_height = exon_height/scaling_factor
                    m_draw_height = e_draw_height
                
                # Finially, draw the damn thing
                # (x, y, width, height)
                # if there are overlaps, the y coordinate gets moved down by a multiple of motif height
                ctx.rectangle(start+edge, m_draw_height+motif_height*prev_overlaps, len(motif_text), motif_height)
                ctx.fill() 
                # after drawing, if this sequence has overlaps, increment prev_overlaps, otherwise reset
                if num_overlaps > 0:
                    # if 
                    prev_overlaps += 1
                else:
                    prev_overlaps = 0

        # export drawing
        surface.write_to_png(self.name+".png")


# parse args
args = parse_args()

# Read in the files 
seqs, names = parse_fasta(args.f)
mots = read_motifs(args.m)

# initalize sequence objects for each sequence, store in list 
seq_objects = []
for seq, name in zip(seqs, names):
    # get starting and ending indices of exon
    match = re.search("[a-z][A-Z]",seq)
    # first capital letter is end of span
    start = match.span()[1]
    match = re.search("[A-Z][a-z]",seq)
    # first lowercase letter is end of span
    end = match.span()[1]
    seq_instance = Sequence(exon_start_=start, exon_stop_=end, motifs_list=mots, sequence=seq, name_=name)
    seq_objects.append(seq_instance)

# Search all of them for all motifs
for seq in seq_objects:
    seq.search()
    #print(seq.found)

# refactor found dict to include whether it overlaps
for seq in seq_objects:
    seq.refactor_found()
    #print(seq.overlap)


# Draw
drawing = Drawing(seq_objects, "fig_", motifs=mots)
drawing.draw_seqs()