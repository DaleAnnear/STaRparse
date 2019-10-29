#!/usr/bin/env python3

def main():
    from TRF_Filter import construct

    if args.output[len(args.output)-1] != "/":
        outpath = args.output+"/"
    else:
        outpath = args.output

    alignment_paramters = [args.match, args.nonmatch, args.gap, args.extendgap]
    for x in alignment_paramters:
        if not x:
            alignment_paramters = [1/2, -2, -5, -1]

    names = ["Chr", "Start", "End", "Pattern_Size", "Units", "Copies_Aligned", "Match%", "Indels%", "Alignment_Score", "A%", "C%", "G%", "T%", "Entropy", "Pattern", "Motif"]
    content = {"Chr":[], "Start":[], "End":[], "Pattern_Size":[], "Units":[], "Copies_Aligned":[], "Match%":[], "Indels%":[], "Alignment_Score":[], "A%":[], "C%":[], "G%":[], "T%":[], "Entropy":[], "Pattern":[], "Motif":[]}

    with open(args.input) as f:
        for line in f:
            for x in range(len(names)):
                content[names[x]].append(line.strip().split()[x])

    repeatdf = construct(content, names, alignment_paramters)
    out_code = "".join([str(datetime.datetime.today().strftime('%Y-%m-%d')), "_", str(random.randint(1,10000))])

    outfile = outpath + "TRF_Panel_Filtered_" + out_code + ".bed"

    repeatdf.to_csv(outfile, index = None, header=False, sep="\t")
    print("File written:   ", outfile)

    path2script = '/home/dannear/STaRparse/bin/Source/Scripts/repeatpanel/TRF_Construct.R'
    arguments = [outpath, outfile, out_code]
    cmd = ['Rscript', path2script] + arguments
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.stdout.read()

    print("File Written:   ", outpath + "TRF_Panel_GangSTR_" + out_code + ".bed")
    print("File Written:   ", outpath + "TRF_Panel_ExpansionHunter_" + out_code + ".json")
    print("###############     COMPLETE     ###############")
