#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:19:37 2019

@author: dannear
"""

import pandas as pd, shlex, subprocess

def annovar(annodf, out_code):
    print("###############     ANNOTATING CGG-REPEAT LOCATIONS    ###############")
    annodf.to_csv('/home/dannear/STaRparse/Annovar/Region_Inputs/Annovar_Input_'+out_code+'.csv', sep = '\t', index = False, header = False)
    args_str = "perl /home/dannear/annovar/annotate_variation.pl -out /home/dannear/STaRparse/Annovar/Output/Annovar_Output_"+out_code+" -build hg19 /home/dannear/STaRparse/Annovar/Region_Inputs/Annovar_Input_"+out_code+".csv /home/dannear/annovar/humandb"
    args = shlex.split(args_str)
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    print(p.stdout.read())
    genes = pd.read_csv("/home/dannear/STaRparse/Annovar/Output/Annovar_Output_"+out_code+".variant_function", sep = '\t', header = None)
    return genes