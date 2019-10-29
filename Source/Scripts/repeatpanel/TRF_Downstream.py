#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:36:57 2019

@author: dannear
"""

import pandas as pd
from math import ceil
from Bio import Align, SeqIO, pairwise2
from Bio.Seq import Seq

names = ["Chr", "Start", "End", "Pattern_Size", "Units", "Copies_Aligned", "Match%", "Indels%", "Alignment_Score", "A%", "C%", "G%", "T%", "Entropy", "Pattern", "Motif"]
content = {"Chr":[], "Start":[], "End":[], "Pattern_Size":[], "Units":[], "Copies_Aligned":[], "Match%":[], "Indels%":[], "Alignment_Score":[], "A%":[], "C%":[], "G%":[], "T%":[], "Entropy":[], "Pattern":[], "Motif":[]}

with open("/media/dannear/Storage/TRF_Data/CGG_Repeat_Panel_04092019.bed")as f:
    for line in f:
        for x in range(len(names)):
            content[names[x]].append(line.strip().split()[x])

def HammingDistance(String1, String2):
    count = 0
    for i in range(0, len(String1)):
       if String1[i] != String2[i]:
           count += 1
    return count

def split_pattern(pattern, motif):
    rep = 0
    indel = 0
    new_mots = {"Mot":[], "Pos":[]}
    start = 0
    counter = False

    for y in range(0, len(motif), 3):

        if indel >= rep*1.5 and rep != 0:
            if len(motif[start:y])/3 > 4:
                new_mots["Pos"].append(y-4)
                new_mots["Mot"].append(motif[start:y-3])
            rep = 0
            counter = True
            indel = 0

        ham = HammingDistance(motif[y:y+3], pattern)
        if ham == 0:
            match = True
            rep += 1
        else:
            match = False
            if ham == 1: indel += 2
            elif ham == 2: indel += 4
            elif ham == 3: indel += 6

        if match == True and counter == True:
            rep = 0
            indel = 0
            start = y
            counter = False

    return new_mots


def alignment_split(pattern, motif, ori):
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

        alignments = pairwise2.align.localms(seq1, seq2, 1/2, -2, -5, -1)
        index = len(alignments)-1
        if index > -1:
            start = alignments[index][3]                # ? conditional required
            end = alignments[index][4]                  # ? conditional rquired
            spliter = motif[start:end]
            if len(spliter) > 10:
                new_motifs["Mot"].append(spliter)
                new_motifs["Pos"].append(int(ori)+int(start)+1)
            sequences.extend(motif.split(spliter, 2))
        else:
            break

        return new_motifs



new_reads = {"Chr":[], "Start":[], "End":[], "Pattern_Size":[], "Units":[], "Copies_Aligned":[], "Match%":[], "Indels%":[], "Alignment_Score":[], "A%":[], "C%":[], "G%":[], "T%":[], "Entropy":[], "Pattern":[], "Motif":[]}
to_pop = []
for x in range(len(content["Motif"])):
    if float(content["Match%"][x]) < 80:
        new_motifs = alignment_split(content["Pattern"][x], content["Motif"][x], content["Start"][x])
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
final.to_csv(r'/media/dannear/Storage/tmp/test_new_bed.bed', index = None, header=False, sep="\t")
l_df.to_csv(r'/media/dannear/Storage/tmp/left_overs.bed', index = None, header=False, sep="\t")
#GGCGGCGACGGCGGTGGCGGCGT GGC GGC GAC GGC GGT GGC GGC GTC