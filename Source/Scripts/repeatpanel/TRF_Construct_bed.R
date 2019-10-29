#!/usr/bin/env Rscript

library(dplyr)

bedfile <- read.csv("/home/dannear/STaRparse/bin/Source/Scripts/genepanel/TRF_CGG_Repeat_Panel.bed", sep = "\t", header = FALSE)

gang_bed <- select(bedfile, c("V1", "V2", "V3", "V4", "V15"))

write.table(gang_bed, "/media/dannear/Storage/tmp/TRF_panel_gangSTR.bed", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


