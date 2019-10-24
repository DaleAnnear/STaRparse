#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

#       Imported merged CSV file
csv <- read.csv("/media/dannear/Storage/csvmerge_test.csv", sep = "\t")

#       Import repeat bed file
ref_bed <- read.csv("/media/dannear/Storage/TRF_Data/CGG_Repeat_Panel_04092019.bed", header = F, sep = "\t")
nam <- c("Chr", "Start", "End", "Pattern_Size", "Units", "Copies_Aligned", "Match%", "Indels%", "Alignment_Score", "A%", "C%", "G%", "T%", "Entropy", "Pattern", "Motif")
names(ref_bed) <- nam
ref_bed$Call_ID <- as.factor(gsub("chr", "", as.character(paste(ref_bed$Chr, ref_bed$Start, sep = "."))))

#       Generate consensus allele lengths
#       Collect all equal reads
ran <- 5
basedf <- csv[ which( between(csv$Allele1.Gang, csv["Allele1.ExH"]-ran, csv["Allele1.ExH"]+ran) |  
                       between(csv$Allele1.Gang, csv["Allele2.ExH"]-ran, csv["Allele2.ExH"]+ran) | 
                        between(csv$Allele2.Gang, csv["Allele1.ExH"]-ran, csv["Allele1.ExH"]+ran) | 
                          between(csv$Allele2.Gang, csv["Allele2.ExH"]-ran, csv["Allele2.ExH"]+ran) ), ]
nadf <- rbind(csv[which(is.na(csv$Ref.ExH)),], csv[which(is.na(csv$Ref.Gang)),])
nadf$Allele1.ExH[which(is.na(nadf$Allele1.ExH))] <- nadf$Allele1.Gang[which(is.na(nadf$Allele1.ExH))]
nadf$Allele2.ExH[which(is.na(nadf$Allele2.ExH))] <- nadf$Allele2.Gang[which(is.na(nadf$Allele2.ExH))]
nadf$Ref.ExH[which(is.na(nadf$Ref.ExH))] <- nadf$Ref.Gang[which(is.na(nadf$Ref.ExH))]

mergedf <- rbind(basedf, nadf)

#       Select out differing reads
x <- rbind(csv, mergedf)
x <- x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(csv), ]
x <- merge(x, ref_bed[c("Call_ID", "Match%", "Entropy")])
x <- na.omit(x)

#       Select out reference pure reads
match <- x[which(x$`Match%` >= 85),]
y <- rbind(x, match)
y <- y[! duplicated(y, fromLast=TRUE) & seq(nrow(y)) <= nrow(x), ]

#       Select out low unit repeats
lowend <- y[which(y$Ref.ExH < 9),]
z <- rbind(y, lowend)
z <- z[! duplicated(z, fromLast=TRUE) & seq(nrow(z)) <= nrow(y), ]

#       Build consensus list
compiledf <- rbind(mergedf, select(match, -c("Match%", "Entropy")), select(lowend, -c("Match%", "Entropy")))
compiledf <- select(compiledf, -c("Ref.Gang", "Allele1.Gang", "Allele2.Gang"))
names(compiledf) <- c('Call_ID', 'Sample_ID', 'Chr', 'Start', 'Ref_Units', 'Allele1_Units', 'Allele2_Units')
z <- select(z[which(z$Allele1.Gang !=1 | z$Allele2.Gang != 1),], -c("Ref.Gang", "Allele1.ExH", "Allele2.ExH", "Match%", "Entropy"))
names(z) <- c('Call_ID', 'Sample_ID', 'Chr', 'Start', 'Ref_Units', 'Allele1_Units', 'Allele2_Units')

consensusdf <- rbind(compiledf, z)
write.table(consensusdf, file="/media/dannear/Storage/tmp/test_consensus.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
View(consensusdf)
