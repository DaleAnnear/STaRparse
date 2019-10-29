#!/usr/bin/env python3

def main():
    path2script = '/home/dannear/STaRparse/bin/Source/Scripts/summaries/Summaries.R'
    arguments = [args.csvinput, args.summaryoutput]
    cmd = ['Rscript', path2script] + arguments
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.stdout.read()

    print("Summary File Written, Sample        : " + args.summaryoutput + "_by_Sample.csv")
    print("Summary File Written, Chromosome    : " + args.summaryoutput + "_by_Chr.csv")
    print("Summary File Written, Locus         : " + args.summaryoutput + "_by_Locus.csv")
    print("Summary File Written, Stability     : " + args.summaryoutput + "_by_Stability.csv")
    print("###############     COMPLETE     ###############")