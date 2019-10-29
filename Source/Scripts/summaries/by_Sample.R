#!/usr/bin/env Rscript

by_sample = function(df){
  medi <- aggregate(list("Median"=c(df$All1,df$All2)), by = list("CallID"=c(df$CallID,df$CallID)), FUN=median)
  df <- merge(df, medi, by="CallID")
  poly <- subset(df, All1 != Median | All2 != Median)
  pa <- aggregate(poly$SampleID, by = list(poly$SampleID), FUN = length)
  a <- aggregate(df$SampleID, by = list(df$SampleID), FUN = length)
  b <- aggregate((All1+All2)/2 ~SampleID, data = df, FUN = mean)
  c <- aggregate(max_reps$MaxReps, by = list(max_reps$SampleID), FUN = max)
  z <- data.frame("SampleID"=a$Group.1, "Highest_Repeat"=c$x, "Total_Repeats"=a$x, "Mean_Rep_length"=round(b$"(All1 + All2)/2", 2), "Unstable_Reps"=pa$x, "Percentage_Unstable_Reps"=round((pa$x/a$x*100), 2))
  z <- z[order(z$SampleID),]
  return(z)}

myArgs <- commandArgs(trailingOnly = TRUE)
readsdf <- read.csv(myArgs[1], header = TRUE, sep ="\t")
names(readsdf) <- c("Call_ID", "Sample_ID", "Chr", "Start", "Ref_Units", "All1", "All2", "Region", "Gene")
bylocusdf <- by_sample(readsdf)

write.table(bylocusdf, paste0(myArgs[2], "_By_Sample.csv"), col.names = TRUE, sep = "\t", row.names = FALSE)



