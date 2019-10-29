#!/usr/bin/env python3

def to_DB(df, dbname, out_code):
#           EXPORT TO DATABASE AND CSV
    print("###############     EXPORTING TO DATABASE    ###############")
    engine = sqlalchemy.create_engine("mysql+pymysql://dannear:3vVrnhJ3@143.169.238.18/dannear")
    df.to_sql("CGG_Repeats_"+dbname+"_"+out_code, engine, if_exists='replace', index=False)

def Output(vcf_data, out_code, output, db):
    vcf_data.to_csv(output+"CGG_Repeats_"+db+"_"+out_code+".csv", sep=',', index=None, header=True)
    if db != "":
        to_DB(vcf_data, db, out_code)
    else:
        print("WARNING: \t VCF data will not be exported to a database. No database was specified.")
    print("###############     COMPLETE     ###############")
    print("Output CSV file can be found at: "+output+"CGG_Repeats_"+db+"_"+out_code+".csv")

def main():
    if args.output[-1] != "/":
        args.output += "/"
    if args.vcfinput[-1] != "/":
        args.vcfinput += "/"
    out_code = "".join([str(datetime.datetime.today().strftime('%Y-%m-%d')), "_", str(random.randint(1,10000))])
    files = glob.glob(args.vcfinput+'*.vcf')

    if files == []:
        print("WARNING: \t No VCF files were found at the specified directory.")
        exit(1)
    else:
        print(files)
        if args.database is None:
            db = ""
        else:
            db = args.database

    if args.STR_Genotyper == "E":
        out_code = "ExpansionHunter_" + out_code
        from From_ExpansionHunter import ExpansionHunter
        vcf_data = ExpansionHunter(files, out_code)
    if args.STR_Genotyper == "G":
        out_code = "GangSTR_" + out_code
        from From_GangSTR import GangSTR
        vcf_data = GangSTR(files, out_code)
    Output(vcf_data, out_code, args.output, db)
