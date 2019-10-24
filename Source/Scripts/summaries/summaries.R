library(dplyr)

savename = "GangSTR_Anno"
csvpath = "/media/dannear/Storage/GangSTR_Data/Results/CSV_Files/CGG_Repeats_GangSTR_2019-07-16_5417.csv"
gangstrpath <- "/media/dannear/Storage/GangSTR_Data/Results/Data_Analysis/"
col_titles <-c("CallID", "SampleID", "Chr", "Start", "End", "GT", "Ref_Units", "All1", "All2", "Alt_STR", "Ref_STR", "Region", "Gene")
all <- read.csv(csvpath)
all$Gene <- gsub(",","->", fixed=T, gsub("\\([^)]+?\\)", "", all$Gene, perl=T))
names(all) <- col_titles
all$CallID <- as.vector(all$CallID)

Sample_Summary <- by_sample(all)
write.csv(Sample_Summary,paste(gangstrpath,"Sample_Summary - ",savename,".csv", sep = ""))
View(Sample_Summary)

Locus_Summary <- by_locus(all)
write.csv(Locus_Summary,paste(gangstrpath,"Locus_Summary - ",savename,".csv", sep = ""))
View(Locus_Summary)

Chr_Summary <- by_Chr(all, Locus_Summary)
write.csv(Chr_Summary,paste(gangstrpath,"Chr_Summary - ",savename,".csv", sep = ""))
View(Chr_Summary)

rat <- ratio_of_stability(Locus_Summary)
write.csv(rat,paste(gangstrpath, "Ratio_of_Stability - ",savename,".csv", sep = ""))
View(rat)

genes <- genelist(Locus_Summary)
write.table(genes,paste(gangstrpath, "Gene List - ",savename,".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

#     RATIO OF STABLE REPEATS VS UNSTABLE REPEATS BY AVERAGE REPEAT LENGTH
ratio_of_stability = function(df) {
  rat <- df[c("CallID", "Mean_Units", "Med_Units", "Instability_Rating", "Status")]
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


#     AGGREGATE BY SAMPLE
by_sample = function(df) {
  df <- df[c("CallID", "SampleID", "Chr", "Start", "GT", "Ref_Units", "All1", "All2", "Gene")]
  medi <- aggregate(list("Median"=c(df$All1,df$All2)), by = list("CallID"=c(df$CallID,df$CallID)), FUN=median)
  df <- merge(df, medi, by="CallID")
  poly <- subset(df, All1 != Median | All2 != Median)

  pa <- aggregate(poly$SampleID, by = list(poly$SampleID), FUN = length)
  a <- aggregate(df$SampleID, by = list(df$SampleID), FUN = length)
  b <- aggregate((All1+All2)/2 ~SampleID, data = df, FUN = mean)
  c <- aggregate(max_reps$MaxReps, by = list(max_reps$SampleID), FUN = max)
  z <- data.frame("SampleID"=a$Group.1, "Highest_Repeat"=c$x, "Total_Repeats"=a$x, "Mean_Rep_length"=round(b$"(All1 + All2)/2", 2), "Unstable_Reps"=pa$x, "Percentage_Unstable_Reps"=round((pa$x/a$x*100), 2))
  z <- z[order(z$SampleID),]
return(z) }

#     AGGREGATE BY LOCUS
by_locus = function(df){
  #df <- df[c("CallID", "SampleID", "Chr", "Start", "GT", "Ref_Units", "All1", "All2")]
  df$avealleles <- (df$All1+df$All2)/2
  avera <- aggregate(list("Mean"=df$avealleles), by = list("CallID"=df$CallID), FUN = mean)
  avera <- avera[order(avera$CallID),]
  medi <- aggregate(list("Median"=df$avealleles), by = list("CallID"=df$CallID), FUN = median)
  medi <- medi[order(medi$CallID),]
  total <- merge(df,medi,by="CallID")
  total$stab <- ifelse(total$All1!=total$Median & total$All2!=total$Median, 1, 
                         ifelse(total$All1==total$Median & total$All2!=total$Median, 0.5, 
                                ifelse(total$All1!=total$Median & total$All2==total$Median, 0.5, 0)))
  ref <- aggregate(list("Ref_Units"=total$Ref_Units), by = list("CallID"=total$CallID, "Chr"=total$Chr, "Start"=total$Start), FUN = mean)
  counts <- aggregate(list("Count"=df$Call), by = list("CallID"=df$CallID), FUN = length)
  stability <- aggregate(list("Stab"=total$stab), by = list("CallID"=total$CallID), FUN = mean)
  maxreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("CallID"=df$CallID), FUN = max)
  minreps <- aggregate(list("All1"=df$All1, "All2"=df$All2), by = list("CallID"=df$CallID), FUN = min)
  stddev <- aggregate(list("All"=c(df$All1,df$All2)), by = list("CallID"=c(df$CallID,df$CallID)), FUN=sd)
  iqr <- aggregate(list("IQR"=c(df$All1,df$All2)), by = list("CallID"=c(df$CallID,df$CallID)), FUN=IQR)
  gen <- df %>% distinct(CallID, Region, Gene)
  gen$CallID <- as.vector(gen$CallID) 
  gen <- gen[order(gen$CallID),]
  counts <- counts[order(counts$CallID),]
  ref <- ref[order(ref$CallID),]
  stability <- stability[order(stability$CallID),]
  maxreps <- maxreps[order(maxreps$CallID),]
  minreps <- minreps[order(minreps$CallID),]
  stddev <- stddev[order(minreps$CallID),]
  iqr <- iqr[order(iqr$CallID),]
  compile <- data.frame("CallID"=ref$CallID, "Chr"=ref$Chr, "Start"=ref$Start, "Ref_Units"=ref$Ref_Units, "Mean_Units"=round(avera$Mean), "Med_Units"=round(medi$Median), "Min_Rep"=pmin(minreps$All1, minreps$All2), "Max_Rep"=pmax(maxreps$All1, maxreps$All2), "Hits"=counts$Count, "Region"=gen$Region, "Gene"=gen$Gene, "Standard_Dev"=stddev$All,"Instability_Rating"=stability$Stab, "IQR"=iqr$IQR)
  compile$Status <- ifelse((compile$Instability_Rating == 0 & compile$Ref_Units != compile$Mean_Units), "FLAGGED", ifelse(compile$Instability_Rating == 0, "STABLE", "UNSTABLE"))
  ev <- expansion_event(df, compile)
  compile$Expansion_Event <- ev$Expansion_Event
  return(compile)} 

#     AGGREGATE BY CHROMOSOME
by_Chr = function(df, locus){
  chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  mbp <- c(249,237,192,183,174,165,153,135,132,132,132,123,108,105,99,84,81,75,69,63,54,57,141,60) 
  chrmbp = data.frame("Chr"=chr, "Mbp"=mbp)
  tot <- aggregate(list("Reps"=locus$Chr), by = list("Chr"=locus$Chr), FUN = length)
  tot$Chr <- ifelse(tot$Chr=="X", 23, ifelse(tot$Chr=="Y", 24, tot$Chr))
  tot <- tot[order(tot$Chr),]
  ave <- aggregate(list("Ave"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = mean)
  ave$Chr <- ifelse(ave$Chr=="X", 23, ifelse(ave$Chr=="Y", 24, ave$Chr))
  ave <- ave[order(ave$Chr),]
  med <- aggregate(list("Med"=(df$All1+df$All2)/2), by = list("Chr"=df$Chr), FUN = median)
  med$Chr <- ifelse(med$Chr=="X", 23, ifelse(med$Chr=="Y", 24, med$Chr))
  med <- med[order(med$Chr),]
  status <- subset(locus, Status == "STABLE")
  stab <- aggregate(list("Status"=status$Status), by = list("Chr"=status$Chr), FUN = length)
  status <- subset(locus, Status != "STABLE")
  instab <- aggregate(list("Status"=status$Status), by = list("Chr"=status$Chr), FUN = length)
  stab <- stab[order(stab$Chr),]
  #instab[nrow(instab) + 1,] = list("Y",0)
  instab <- instab[order(instab$Chr),]
  z <- data.frame("Chromosome"=chrmbp$Chr, "Number_of_Repeats"=tot$Rep, "Mean_Rep_Len"=ave$Ave, "Median_Rep_Len"=med$Med, "Repeats_per_Mbp"=round(tot$Reps/chrmbp$Mbp, 3), "Poly_per_Mbp"=round(instab$Status/chrmbp$Mbp, 3), "Stable_Reps"=stab$Status, "Unstable_Reps"=instab$Status)
return(z)}

#     Significant Expansion Event
expansion_event = function(alldf, loci){
  loci$CI_low <- loci$Mean_Units - (qnorm(0.95)*(loci$Standard_Dev))-1
  loci$CI_upp <- loci$Mean_Units + (qnorm(0.95)*(loci$Standard_Dev))+1
  loci$Expansion_Event <- ifelse(loci$Min_Rep < loci$CI_low & loci$Max_Rep > loci$CI_upp, "EC", ifelse(loci$Max_Rep > loci$CI_upp, "E", ifelse(loci$Min_Rep < loci$CI_low, "C", "S")))
  out <- loci[c("CallID", "Expansion_Event")]
return(out)}

#     Get Gene List
genelist = function(df){
  genes <- gsub(",","->", fixed=T, gsub("\\([^)]+?\\)", "", df$Gene,perl=T))
  genes <- unique(unlist(strsplit(genes, "->")))
return(genes)}

