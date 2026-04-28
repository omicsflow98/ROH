#!/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

option_list = list(
  make_option(c("-T", "--tsv"), default=NULL,
              help = "merged tsv file", metavar = "tsv file"),
  make_option(c("-V", "--version"), default=NULL,
              help = "version of the nextflow process", metavar = "process version"),
  make_option(c("-L", "--length"), type = "numeric", default=NULL,
              help = "Genome length in MB", metavar = "Genome Length"),
  make_option(c("-S", "--state"), default=NULL,
              help = "HOM or HET", metavar = "analysis type")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

if(is.null(opt$tsv)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call. = FALSE)

}

file <- read.table(opt$tsv, sep = "\t", fill = TRUE, header = TRUE) %>%
  arrange(CHR, POS1, POS2)

indiv_avg <- file %>%
  mutate(IID = as.factor(IID)) %>%
  select(IID, KB) %>%
  rename(KBold = KB) %>%
  group_by(IID) %>%
  summarise(NSEG = n(),
    KB = round(sum(KBold, na.rm = TRUE), 1),
    KBAVG = round(mean(KBold), 2)) %>%
  mutate(FRR = round(KB / (opt$length*1000), 3)) %>%
  ungroup()

 write.table(indiv_avg,file=paste0("indiv.tsv"),quote = FALSE,row.names = FALSE,sep = "\t")

result_temp <- file %>%
  group_by(CHR, POS1, POS2) %>%
  summarise(
    samples = paste(IID, collapse = ","),
    CHR = first(CHR),
    across(KB:last_col(), ~first(.)),
    .groups = "drop"
  )

  if (opt$state == "HOM") {
    result <- result_temp %>%
      mutate(ROH_ID = paste0("ROH_", row_number())) %>%
      select(ROH_ID, POS1, POS2, samples, everything())
  } else {
    result <- result_temp %>%
      mutate(HRR_ID = paste0("HRR_", row_number())) %>%
      select(HRR_ID, POS1, POS2, samples, everything())
  }


i <- 0

for (chrom in unique(result$CHR)) {
  i <- i + 1
  one_chrom <- result %>%
    filter(CHR == chrom) %>%
    arrange(POS1, POS2)

  write.table(one_chrom,file=paste0("file_", i, ".tsv"),quote = FALSE,row.names = FALSE,sep = "\t")
}

ver <- strsplit(R.version.string, " ")[[1]][3]
settings_list <- list(
  `Process version` = list(
    splitROH = opt$version
  ),
  `Tool version` = list(
    R =  ver
  )
)

write_yaml(settings_list, file = "versions.yml")
