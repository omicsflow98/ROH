#!/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-T", "--tsv"), default=NULL,
              help = "merged tsv file", metavar = "tsv file")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

if(is.null(opt$tsv)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call. = FALSE)

}

file <- read.table(opt$tsv, sep = "\t", fill = TRUE, header = TRUE) %>%
  arrange(CHR, POS1, POS2)

result <- file %>%
  group_by(CHR, POS1, POS2) %>%
  summarise(
    samples = paste(FID, collapse = ","),
    CHR = first(CHR),
    across(KB:last_col(), ~first(.)),
    .groups = "drop"
  ) %>%
  mutate(ROH_ID = paste0("ROH_", row_number())) %>%
  select(ROH_ID, POS1, POS2, samples, everything())

i <- 0

for (chrom in unique(result$CHR)) {
  i <- i + 1
  one_chrom <- result %>%
    filter(CHR == chrom) %>%
    arrange(POS1, POS2)

  write.table(one_chrom,file=paste0("file_", i, ".tsv"),quote = FALSE,row.names = FALSE,sep = "\t")
}