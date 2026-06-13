#!/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(plotly))
suppressMessages(library(DT))

subpopulations <- read.table(subpop, 
                             sep = "\t", 
                             fill = TRUE, 
                             header = TRUE,
                             check.names = FALSE)
colnames(subpopulations) <- gsub("\\s+", "_", colnames(subpopulations))

#number of samples per subpopulation
num_samples <- colSums(subpopulations[,-1])
num_samples_df <- enframe(num_samples) %>%
  rename(Subpopulation = name, `Number of Samples` = value)
num_samples_df$Subpopulation <- gsub("_", " ", num_samples_df$Subpopulation)

index <- which(subpopulations == 1, arr.ind = TRUE) 
names <- colnames(subpopulations)[index[, "col"]]
subpopulations[index] <- names
subpopulations[subpopulations == 0] <- NA 
subpopulations <- subpopulations %>%
  unite("groupings", -SampleID, na.rm = TRUE, sep = ",", remove = FALSE) %>%
  select(SampleID, groupings)

if (hom) {
  ROHfile <- read.table(indiv_homfile, sep = "\t", fill = TRUE, header = TRUE)
  
  merged_ROH <- merge(ROHfile, subpopulations, by.x = "IID", by.y = "SampleID")

  #number of ROH in each subpopulation
#  plot_nROH <- merged_ROH %>%
#    select(groupings, NSEG) %>%
#    group_by(groupings) %>%
#    summarise(value = sum(NSEG))
#  plot_nROH$groupings <- gsub("_", " ", plot_nROH$groupings)

plot_nROH <- merged_ROH %>%
  select(groupings, NSEG)
  
  #ROH Boxplot
  plot_violin_ROH <- merged_ROH %>%
  select(groupings, KBAVG) 
}

if (het) {
  HRRfile <- read.table(indiv_hetfile, sep = "\t", fill = TRUE, header = TRUE)
  
  merged_HRR <- merge(HRRfile, subpopulations, by.x = "IID", by.y = "SampleID")

  #number of ROH in each subpopulation
#  plot_nHRR <- merged_HRR %>%
#    select(groupings, NSEG) %>%
#    group_by(groupings) %>%
#    summarise(value = sum(NSEG))
#  plot_nHRR$groupings <- gsub("_", " ", plot_nHRR$groupings)

plot_nHRR <- merged_HRR %>%
  select(groupings, NSEG)
  
  #ROH Boxplot
  plot_violin_HRR <- merged_HRR %>%
  select(groupings, KBAVG) 
}
