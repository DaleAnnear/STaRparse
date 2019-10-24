#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:21:44 2019

@author: dannear
"""
import pandas as pd, os, vcf, sys

def GangSTR(vcffiles, out_code):
#           IMPORT VCF FILES AND DESIRED EXTRACT DATA
    print("###############     IMPORTING VCF DATA     ###############")
    mdf = {'Call':[], 'Sample_ID':[], 'Chr':[], 'Start':[], 'End':[], 'GT':[], 'Ref_Units':[], "Allele1_Units":[], "Allele2_Units":[], "Alternate_STR":[], "Reference_STR":[]}
    anno = {'Chr':[], 'PoS':[], 'PoE':[], 'RA':[], 'AA':[]}
    count = 0
    for x in vcffiles:
        count += 1
        base = os.path.splitext(os.path.basename(x))[0]
        print("Processing VCF file: ", base, "\tFile no.: ", count)
        my_vcf = vcf.Reader(filename=x)
        for record in my_vcf:
            samp = record.samples
            for i in samp:
                 if i['GT'] == ".":
                     continue
                 if record.CHROM[3:] in ["X", "Y"]:
                     if record.CHROM[3:] == "X": Call = "X."+str(record.POS)
                     if record.CHROM[3:] == "Y": Call = "Y."+str(record.POS)
                 else:
                     Call = str(record.CHROM[3:])+"."+str(record.POS)
                 chrom = record.CHROM[3:]
                 start = int(record.POS)
                 end = int(record.INFO["END"])
                 if i['GT'] in ["0/1", "1/0", "1/1", "1/2"]:
                     units = i["REPCN"]
                 else:
                     units = [round(record.INFO["REF"]),round(record.INFO["REF"])]
                 mdf['Call'].append(Call)
                 mdf['Sample_ID'].append(str(base.split('.bam')[0]))
                 mdf['Chr'].append(chrom)
                 mdf['Start'].append(start)
                 mdf['End'].append(end)
                 mdf['GT'].append(i['GT'])
                 mdf['Ref_Units'].append(round(record.INFO["REF"]))
                 mdf["Allele1_Units"].append(units[0])
                 mdf["Allele2_Units"].append(units[1])
                 mdf["Alternate_STR"].append(str(record.ALT)[1:-1])
                 mdf["Reference_STR"].append(record.REF)

#           ANNOTATE LOCI WITH ANNOVAR
                 anno['Chr'].append(chrom)
                 anno['PoS'].append(start)
                 anno['PoE'].append(end)
                 anno['RA'].append(0)
                 anno['AA'].append("-")

    annodf = pd.DataFrame(anno)
    from Run_ANNOVAR import annovar
    genes = annovar(annodf, out_code)
    df = pd.DataFrame(mdf)
    df['Region'] = pd.Series(genes[0])
    df['Gene'] = pd.Series(genes[1])
    return df