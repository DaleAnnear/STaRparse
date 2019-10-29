#!/usr/bin/env Rscript

suppressMessages(library(dplyr))

myArgs <- commandArgs(trailingOnly = TRUE)

#         LOAD FUNCTIONs
#         JSON
con_json = function(outpath, infile, outcode) {
  input <- read.csv(infile, sep = "\t", header = FALSE) #, col.names = c("chr", "start", "end", "X3", "motif"))
  bed <- select(input, c("V1", "V2", "V3", "V15"))
  names(bed) <- c("chr", "start", "end", "motif")
  bed$chr <- gsub("chr", "", bed$chr)
  bed$JSON <- sprintf(
    "{
    \"LocusId\": \"%s\", 
    \"LocusStructure\": \"(%s)*\", 
    \"ReferenceRegion\": \"%s:%s-%s\",
    \"VariantType\": \"RareRepeat\"
  },"
    , paste(bed$chr, ".", as.character(bed$start), sep=""), bed$motif, bed$chr, bed$start, bed$end)
  ind <- length(bed$JSON)
  bed$JSON[ind] <- substr(bed$JSON[ind], 1, nchar(bed$JSON[ind])-1)
  bed$JSON[ind] <- paste0(bed$JSON[ind], "]")
  bed$JSON[1] <- paste0("[", bed$JSON[1])
  write(bed$JSON, file = paste0(outpath, "TRF_Panel_ExpansionHunter_", outcode, ".json"))
}
#         BED
con_bed= function(outpath, infile, outcode) {
  bedfile <- read.csv(infile, sep = "\t", header = FALSE)
  gang_bed <- select(bedfile, c("V1", "V2", "V3", "V4", "V15"))
  write.table(gang_bed, paste0(outpath, "TRF_Panel_GangSTR_", outcode, ".bed"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  paste0("File written:", "\t", outpath, "TRF_Panel_GangSTR_", outcode, ".bed")
}

#         Define arguments
outpath <- myArgs[1]
infile <- myArgs[2]
outcode <- myArgs[3]

#         Run functions
con_json(outpath, infile, outcode)
con_bed(outpath, infile, outcode)