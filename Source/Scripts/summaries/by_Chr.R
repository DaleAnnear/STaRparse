#!/usr/bin/env Rscript

library(reshape2)
library(dplyr)

by_locus = function(df){
  df$avealleles <- (df$All1+df$All2)/2
  avera <- aggregate(list("Mean"=df$avealleles), by = list("Call_ID"=df$Call_ID), FUN = mean)
  avera <- avera[order(avera$Call_ID),]
  mdata <- melt(df[c("Call_ID", "All2", "All2")], id="Call_ID")
  medi <- aggregate(list("Median" = mdata$value), by=list("Call_ID" = mdata$Call_ID), median)
  medi <- medi[order(medi$Call_ID),]
  total <- merge(df,medi,by="Call_ID")
  total$stab <- ifelse(total$All1!=total$Median & total$All2!=total$Median, 1, 
                       ifelse(total$All1==total$Median & total$All2!=total$Median, 0.5, 
                              ifelse(total$All1!=total$Median & total$All2==total$Median, 0.5, 0)))
  ref <- aggregate(list("Ref_Units"=total$Ref_Units), by = list("Call_ID"=total$Call_ID, "Chr"=total$Chr, "Start"=total$Start), FUN = mean)
  counts <- aggregate(list("Count"=df$Call_ID), by = list("Call_ID"=df$Call_ID), FUN = length)
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

by_sample = function(df, locus){
  #medi <- locus[c("Call_ID")]#aggregate(list("Median"=c(df$All1,df$All2)), by = list("Call_ID"=c(df$Call_ID,df$Call_ID)), FUN=median)
  df <- merge(df, locus[c("Cal_ID", "Med_Units")], by="Call_ID")
  poly <- subset(df, df$All1 != Median | df$All2 != Median)
  View(poly)
  pa <- aggregate(poly$Sample_ID, by = list(poly$Sample_ID), FUN = length)
  a <- aggregate(df$Sample_ID, by = list(df$Sample_ID), FUN = length)
  b <- aggregate((All1+All2)/2 ~Sample_ID, data = df, FUN = mean)
  c <- aggregate(max_reps$MaxReps, by = list(max_reps$Sample_ID), FUN = max)
  z <- data.frame("Sample_ID"=a$Group.1, "Highest_Repeat"=c$x, "Total_Repeats"=a$x, "Mean_Rep_length"=round(b$"(All1 + All2)/2", 2), "Unstable_Reps"=pa$x, "Percentage_Unstable_Reps"=round((pa$x/a$x*100), 2))
  z <- 0 #z[order(z$Sample_ID),]
  return(z)}

readsdf <- read.csv("/media/dannear/Storage/today_consensus.csv", header = TRUE, sep ="\t")
names(readsdf) <- c("Call_ID", "Sample_ID", "Chr", "Start", "Ref_Units", "All1", "All2", "Region", "Gene")

locus <- by_locus(readsdf)

bysam <- by_sample(readsdf)

View(bysam)


