#!/usr/bin/env Rscript

#         LOAD FUNCTION
con_json = function() {
  bed <- read.csv(opt$bed_file, sep = "\t", col.names = c("chr", "start", "end", "X3", "motif"))
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
  write(bed$JSON, file = paste0(opt$json, ".json"))
  sprintf("SUCCESS:   %s was convereted to a JSON file.", opt$bed_file)
}

#         IMPORT LIBRARIES
library(optparse)

#         INITATE ARGUMENTS
option_list <- list(
  make_option(c("-b", "--bed_file"),type = "character", default="/media/dannear/Storage/TRF_Data/GangSTR_22082019.bed",
              help="Provide bed file with the STR locations"),
  make_option(c("-j", "--json"), type="character", default="/media/dannear/Storage/ExpansionHunter_Data/JSON_Files/default_json",
              help="Provide the file path and name for the JSON output file")
)
parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
file <- args$args

#         RUN SCRIPT
con_json()