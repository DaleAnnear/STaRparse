#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Created on Mon Jul 15 14:54:46 2019 @author: dannear"""

#           IMPORT MODULES
import pandas as pd, argparse, glob, os, vcf, datetime, sqlalchemy, shlex, subprocess, random

#           DEFINE MAIN USING ARGPARSE
def main():
    parser = argparse.ArgumentParser(description='Select VCF files to process')
    parser.add_argument('STR_Genotyper', choices=['E','G','L'], metavar="", help='Special testing value')
    parser.add_argument('-i', '--vcfinput', type=str, metavar="", required=True, help='VCF Input File Path')
    parser.add_argument('-o', '--output', type=str, metavar="", nargs="?", const="/home/dannear/STaRparse/Default_Output", required=False, help='CSV Output File Path')
    parser.add_argument('-db', '--database', type=str, metavar="", required=False, help='Enter Database Name')
    args = parser.parse_args()

    if args.output[-1] != "/":
       args.output += "/"
    if args.vcfinput[-1] != "/":
       args.vcfinput += "/"
    out_code = "".join([str(datetime.datetime.today().strftime('%Y-%m-%d')), "_", str(random.randint(1,10000))])
    files = glob.glob(args.vcfinput+'*.vcf')

    if files == []:
        print("WARNING: \t No VCF files were found at the specified directory.")
    else:
        print(files)
        if args.database is None:
            db = ""
        else:
            db = args.database + "_"

        if args.STR_Genotyper == "E":
            out_code = "ExpansionHunter_" + out_code
            vcf_data = ExpansionHunter(files, out_code)
        if args.STR_Genotyper == "G":
            out_code = "GangSTR_" + out_code
            vcf_data = GangSTR(files, out_code)
        Output(vcf_data, out_code, args.output, db)

def Output(vcf_data, out_code, output, db):
    vcf_data.to_csv(output+"CGG_Repeats_"+db+out_code+".csv", sep=',', index=None, header=True)
    if db != "":
        to_DB(vcf_data, db, out_code)
    else:
        print("WARNING: \t VCF data will not be exported to a database. No database was specified.")
    print("###############     COMPLETE     ###############")
    print("Output CSV file can be found at \t "+output+"CGG_Repeats_"+db+out_code+".csv")

#           DEFINE VCF_to_CSV_to_DB SCRIPT
def ExpansionHunter(vcffiles, out_code):
#           IMPORT VCF FILES AND DESIRED EXTRACT DATA
    print("###############     IMPORTING VCF DATA     ###############")
    mdf = {'Call':[], 'Sample_ID':[], 'Chr':[], 'Start':[], 'End':[], 'GT':[], 'Ref_Units':[], "Allele1_Units":[], "Allele2_Units":[]}
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
                 Call = record.INFO['REPID']
                 chrom = record.CHROM[3:]
                 start = int(record.POS)
                 end = int(record.INFO["END"])
                 units = i['REPCN']
                 mdf['Call'].append(Call)
                 mdf['Sample_ID'].append(base.strip('.bam_ExpansionHunter.vcf'))
                 mdf['Chr'].append(chrom)
                 mdf['Start'].append(start)
                 mdf['End'].append(end)
                 mdf['GT'].append(i['GT'])
                 mdf['Ref_Units'].append(round(record.INFO["REF"]))
                 if chrom == "Y":
                     mdf["Allele1_Units"].append(units)
                     mdf["Allele2_Units"].append(0)
                 else:
                     mdf["Allele1_Units"].append(units.split("/")[0])
                     mdf["Allele2_Units"].append(units.split("/")[1])

#           ANNOTATE LOCI WITH ANNOVAR
                 anno['Chr'].append(chrom)
                 anno['PoS'].append(start)
                 anno['PoE'].append(end)
                 anno['RA'].append(0)
                 anno['AA'].append("-")

    annodf = pd.DataFrame(anno)
    genes = annovar(annodf, out_code)
    df = pd.DataFrame(mdf)
    df['Region'] = pd.Series(genes[0])
    df['Gene'] = pd.Series(genes[1])
    return df

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
                 mdf['Sample_ID'].append(base.strip('.bam_GangSTR.vcf'))
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
    genes = annovar(annodf, out_code)
    df = pd.DataFrame(mdf)
    df['Region'] = pd.Series(genes[0])
    df['Gene'] = pd.Series(genes[1])
    return df

def annovar(annodf, out_code):
    print("###############     ANNOTATING CGG-REPEAT LOCATIONS    ###############")
    annodf.to_csv('/home/dannear/STaRparse/Annovar/Region_Inputs/Annovar_Input_'+out_code+'.csv', sep = '\t', index = False, header = False)
    args_str = "perl /home/dannear/annovar/annotate_variation.pl -out /home/dannear/STaRparse/Annovar/Output/Annovar_Output_"+out_code+" -build hg19 /home/dannear/STaRparse/Annovar/Region_Inputs/Annovar_Input_"+out_code+".csv /home/dannear/annovar/humandb"
    args = shlex.split(args_str)
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    print(p.stdout.read())
    genes = pd.read_csv("/home/dannear/STaRparse/Annovar/Output/Annovar_Output_"+out_code+".variant_function", sep = '\t', header = None)
    return genes

def to_DB(df, dbname, out_code):
#           EXPORT TO DATABASE AND CSV
    print("###############     EXPORTING TO DATABASE    ###############")
    engine = sqlalchemy.create_engine("mysql+pymysql://dannear:3vVrnhJ3@143.169.238.18/dannear")
    df.to_sql("CGG_Repeats_"+dbname+"_"+out_code, engine, if_exists='replace', index=False)


#           EXECUTE MAIN
if __name__ == '__main__':
    main()
