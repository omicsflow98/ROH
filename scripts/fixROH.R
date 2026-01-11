#!/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

CombineROH <- function(list1, list2) {
  list1 <- str_split(list1, ",")
  list1 <- unlist(list1)
  list2 <- str_split(list2, ",")
  list2 <- unlist(list2)

  combine <- unique(append(list1, list2))

  final_row <- paste(combine, collapse = ",")

  return(final_row)
}

option_list = list(
  make_option(c("-T", "--tsv"), default=NULL,
              help = "merged tsv file", metavar = "tsv file"),
  make_option(c("-V", "--version"), default=NULL,
              help = "nextflow process version", metavar = "process version")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

if(is.null(opt$tsv)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call. = FALSE)

}

one_chrom <- read.table(opt$tsv, sep = "\t", fill = TRUE, header = TRUE)

if (nrow(one_chrom) < 2) {

  ver <- strsplit(R.version.string, " ")[[1]][3]
  settings_list <- list(
    `Process Version` = list(
      fixROH = opt$version
    ),
    `Tool Version` = list(
      R =  ver
    )
  )

  write_yaml(settings_list, file = "versions.yml")

  filename <- stringr::str_remove(opt$tsv, "\\..*")
  write.table(one_chrom,file=paste0(filename, "_fixed.tsv"),quote = FALSE,row.names = FALSE,sep = "\t")
  quit(save = "no", status = 0)
  
}

for (row in 1:(nrow(one_chrom)-1)) {
  POS1 <- one_chrom$POS1[row]
  POS2 <- one_chrom$POS2[row]
  for (scroll in (row+1):nrow(one_chrom)) {
    ROHsamples <- one_chrom$samples[row]
    ROHadd <- one_chrom$samples[scroll]

    final_row <- CombineROH(ROHsamples, ROHadd)

    POSstart <- one_chrom$POS1[scroll]
    POSend <- one_chrom$POS2[scroll]

    if (POS2 <= POSstart) {
      break
    }

    if (POS2 >= POSend ) {
      one_chrom$samples[scroll] <- final_row
    } else if (POS1 == POSstart && POS2 < POSend){
      one_chrom$samples[row] <- final_row
    }
  }
}

filename <- stringr::str_remove(opt$tsv, "\\..*")
write.table(one_chrom,file=paste0(filename, "_fixed.tsv"),quote = FALSE,row.names = FALSE,sep = "\t")

ver <- strsplit(R.version.string, " ")[[1]][3]
settings_list <- list(
  `Process version` = list(
    fixROH = opt$version
  ),
  `Tool version` = list(
    R =  ver
  )
)

write_yaml(settings_list, file = "versions.yml")