#!/usr/bin/python

# This script will identify pairs of sequences that are within 3 substitution differences. 

# The input file is expected to contain two columns:
# The first column contains the amino acid CDR3 sequences of the screen and reference sequences
# The second column defines whether the sequence is to be screened ("screen") or if it is a reference sequence ("reference")
# An output file will be created defining all the matches within 3 substitution differences, detailing the number of differences between the reference and screen sequence, and the substring that matched between them.
# Written by Rachael Bashford Rogers - rachael.bashford-rogers@bioch.ox.ac.uk


import math
import sys
from collections import defaultdict
import os
import numpy as np

def Get_sequences(file):
    fh=open(file, "r")
    screen = {}
    reference = {}
    for l in fh:
        l=l.strip().split(',')
        if(l[0]!='screen_seqs'):
            if(l[1]=='screen'):screen[l[0]] = 0
            else:reference[l[0]] = 0
    fh.close()
    return(screen, reference)

def Get_sequence_words(screen,n_words):
    words = {}
    for seq in screen:
        l = len(seq)
        if(l<9):w = [seq[0:int(l/2)], seq[int(l/2):l]]
        else:
            w = [seq[0:int(l/3)],seq[int(l/3):int(l*2/3)], seq[int(l*2/3):l]]
        for i in w:
            if(len(w)>=3):
                if(i in words):words[i] = words[i] + [seq]
                else:words[i] = [seq]
    return(words)

def Screen_sequences(screen, reference,file_out):
    print (len(screen), len(reference))
    # set of reference set of half_sequences
    n_words = 2
    words = Get_sequence_words(screen,n_words)
    print (len(words), "words")
    # per screen, test each half sequence
    match_from,match_to,diff,match_ref_substring = [],[],[],[]
    fail_count = 0
    print_times =np.arange(0, len(reference), 100).tolist()
    done = 0
    fh=open(file_out, "w")
    fh.close()
    out, ind = '#ref_seq\tscreen_seq\tdifference\tscreen_substring\n',0
    for ref in reference:
        done = done+1
        if(done in print_times):
            print(done)
        lref = len(ref)
        for w in words: 
            f = ref.find(w)
            if(f!=-1):
                choices = words[w]
                for c in choices:
                    f2 = c.find(w)
                    start, end = f-f2, f-f2+len(c)
                    ref_comp =  ref[start:end]
                    if(len(ref_comp)==len(c)):
                        dif = [i for i in range(len(ref_comp)) if ref_comp[i]!=c[i]]
                        if(len(dif)<=3):
                            match_from,match_to = match_from+[ref],match_to+[c]
                            diff,match_ref_substring = diff+[len(dif)],match_ref_substring+[ref_comp]
                            if(len(match_from)>10):
                                for i in range(len(match_from)):
                                    out=out+"\t".join([match_from[i],match_to[i], str(diff[i]),match_ref_substring[i]])+"\n"
                                    fh=open(file_out,"a")
                                    fh.write(out)
                                    fh.close()
                                    out = ''
                                match_from,match_to,diff,match_ref_substring = [],[],[],[]
                    elif(start>=0):
                        fail_count = fail_count+1
    print ("fail",fail_count)
    return()

def Check_matches(file_out):
    fh=open(file_out,"r")
    counts = {}
    for l in fh:
        if(l[0]!="#"):
            l=l.strip().split()
            score = int(l[2])
            if(score in counts):
                counts[score] = counts[score]+1
            else:counts[score] = 1
    fh.close()
    print ("score\tcount")
    for i in counts:
        print (i,"\t",counts[i])
    return()



###########
file = "/Users/ssammut/B-cell_immunosurveillance_breast_cancer/output/data/screen_pathogen-sequences_input.early.txt"
file_out = "/Users/ssammut/B-cell_immunosurveillance_breast_cancer/output/data/screen_pathogen-sequences_output.early.txt"


screen, reference = Get_sequences(file)
Screen_sequences(screen, reference,file_out)


Check_matches(file_out)




