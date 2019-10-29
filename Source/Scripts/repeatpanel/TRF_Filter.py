#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:36:57 2019

@author: dannear
"""

import pandas as pd
from math import ceil
from Bio import pairwise2
from Bio.Seq import Seq

def alignment_split(pattern, motif, ori, ali):
    new_motifs = {"Mot":[],"Pos":[]}

    sequences = [motif]

    for x in sequences:
        if len(x) < 12:
            pass
        seq1 = Seq(x)
        ref = ""
        for x in range(ceil(len(x)/3)):
            ref = ref + pattern
        seq2 = Seq(ref)

        alignments = pairwise2.align.localms(seq1, seq2, ali[0], ali[1], ali[2], ali[3])
        index = len(alignments)-1
        if index > -1:
            start = alignments[index][3]
            end = alignments[index][4]
            spliter = motif[start:end]
            if len(spliter) > 10:
                new_motifs["Mot"].append(spliter)
                new_motifs["Pos"].append(int(ori)+int(start)+1)
            sequences.extend(motif.split(spliter, 2))
        else:
            break

        return new_motifs

def construct(content, names, ali):
    new_reads = {"Chr":[], "Start":[], "End":[], "Pattern_Size":[], "Units":[], "Copies_Aligned":[], "Match%":[], "Indels%":[], "Alignment_Score":[], "A%":[], "C%":[], "G%":[], "T%":[], "Entropy":[], "Pattern":[], "Motif":[]}
    to_pop = []
    for x in range(len(content["Motif"])):
        if float(content["Match%"][x]) < 80:
            new_motifs = alignment_split(content["Pattern"][x], content["Motif"][x], content["Start"][x], ali)
            to_pop.append(x)
            for y in range(len(new_motifs["Mot"])):
                for z in range(len(names)):
                    if names[z] == "Start":
                        new_reads[names[z]].append(new_motifs["Pos"][y])
                    elif names[z] == "End":
                        new_reads[names[z]].append(new_motifs["Pos"][y] + len(new_motifs["Mot"][y])-1 )
                    elif names[z] == "Motif":
                        new_reads[names[z]].append(new_motifs["Mot"][y])
                    elif names[z] == "Pattern_Size":
                        new_reads[names[z]].append(round(len(new_motifs["Mot"][y])/3, 1))
                    else:
                        new_reads[names[z]].append(content[names[z]][x])

    c_df = pd.DataFrame(data=content)
    l_df = c_df.loc[to_pop]
    c_df = c_df.drop(to_pop)
    n_df = pd.DataFrame(data=new_reads)

    final = c_df.append(n_df)

    return final
