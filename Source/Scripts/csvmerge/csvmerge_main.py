#!/usr/bin/env python3

def main():
    path2script = '/home/dannear/STaRparse/bin/Source/Scripts/csvmerge/CSV_Merge_Consensus.R'
    arguments = [args.expansionhuntercsv, args.gangstrcsv, args.outputcsv]
    cmd = ['Rscript', path2script] + arguments
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.stdout.read()

    print("Output CSV files can be found at: " + args.outputcsv + "_merged.csv & " + args.outputcsv + "_consensus.csv")

    print("###############     COMPLETE     ###############")