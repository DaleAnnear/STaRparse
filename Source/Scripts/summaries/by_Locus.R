#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(reshape2))


by_locus = function(df){
  df$avealleles <- (df$All1+df$All2)/2
  avera <- aggregate(list("Mean"=df$avealleles), by = list("Call_ID"=df$Call_ID), FUN = mean)
  avera <- avera[order(avera$Call_ID),]
  
  #     Calculate median allele length
  mdata <- melt(df[c("Call_ID", "All2", "All2")], id="Call_ID")
  medi <- aggregate(list("Median" = mdata$value), by=list("Call_ID" = mdata$Call_ID), median)
  medi <- medi[order(medi$Call_ID),]
  total <- merge(df,medi,by="Call_ID")
  total$stab <- ifelse(total$All1!=total$Median & total$All2!=total$Median, 1, 
                       ifelse(total$All1==total$Median & total$All2!=total$Median, 0.5, 
                              ifelse(total$All1!=total$Median & total$All2==total$Median, 0.5, 0)))
  ref <- aggregate(list("Ref_Units"=total$Ref_Units), by = list("Call_ID"=total$Call_ID, "Chr"=total$Chr, "Start"=total$Start), FUN = mean)
  counts <- aggregate(list("Count"=df$Call), by = list("Call_ID"=df$Call_ID), FUN = length)
  stability <- aggregate(list("Stab"=total$stab), by = list("Call_ID"=total$Call_ID), FUN = mean)
  maxreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = max)
  minreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("Call_ID"=df$Call_ID), FUN = min)
  stddev <- aggregate(list("Alls"=c(df$All1,df$All2)), by = list("Call_ID"=c(df$Call_ID,df$Call_ID)), FUN=sd)
  gen <- df %>% distinct(Call_ID, Region, Gene)
  gen$Call_ID <- as.vector(gen$Call_ID) 
  gen <- gen[order(gen$Call_ID),]
  counts <- counts[order(counts$Call_ID),]
  ref <- ref[order(ref$Call_ID),]
  stability <- stability[order(stability$Call_ID),]
  maxreps <- maxreps[order(maxreps$Call_ID),]
  minreps <- minreps[order(minreps$Call_ID),]
  stddev <- stddev[order(minreps$Call_ID),]
  compile <- data.frame("Call_ID"=ref$Call_ID, "Chr"=ref$Chr, "Start"=ref$Start, "Gene"=gen$Gene, "Region"=gen$Region, "Ref_Units"=ref$Ref_Units, "Mean_Units"=round(avera$Mean), "Med_Units"=round(medi$Median), "Min_Rep"=pmin(minreps$All1, minreps$All2), "Max_Rep"=pmax(maxreps$All1, maxreps$All2), "Hits"=counts$Count, "Standard_Dev"=stddev$All,"Instability_Rating"=stability$Stab)
  compile$Status <- ifelse(compile$Instability_Rating == 0, "STABLE", "UNSTABLE")
  return(compile)}

ratio_of_stability = function(df) {
  rat <- df[c("Call_ID", "Mean_Units", "Med_Units", "Instability_Rating", "Status")]
  rat$Med_Units <- round(rat$Med_Units)
  unstable_reps <- subset(rat, rat$Status != "STABLE")
  stable_reps <- subset(rat, rat$Status == "STABLE")
  count_stable <- aggregate(list("Stable_Reps"=stable_reps$Med_Units), by = list("Med_Units"=stable_reps$Med_Units), FUN = length)
  count_unstable <- aggregate(list("Unstable_Reps"=unstable_reps$Mean_Units), by = list("Med_Units"=unstable_reps$Med_Units), FUN = length)
  count_total <- merge(count_stable, count_unstable, all = TRUE)
  count_total[is.na(count_total)] <- 0
  count_total$Total_Reps <- count_total$Stable_Reps+count_total$Unstable_Reps
  count_total$Flagged_Unstable <- (count_total$Unstable_Reps)/count_total$Total_Reps
  return(count_total)}

myArgs <- commandArgs(trailingOnly = TRUE)
readsdf <- read.csv(myArgs[1], header = TRUE, sep ="\t")
names(readsdf) <- c("Call_ID", "Sample_ID", "Chr", "Start", "Ref_Units", "All1", "All2", "Region", "Gene")
bylocusdf <- by_locus(readsdf)
bystab <- ratio_of_stability(bylocusdf)

write.table(bylocusdf, paste0(myArgs[2], "_By_Locus.csv"), col.names = TRUE, sep = "\t", row.names = FALSE)
write.table(bystab, paste0(myArgs[2], "_By_Stability.csv"), col.names = TRUE, sep = "\t", row.names = FALSE)



