#!/usr/bin/env Rscript

suppressMessages(library(detectRUNS))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-P", "--ped"), default=NULL,
              help = "plink PED file", metavar = "PED file"),
  make_option(c("-M", "--map"), default=NULL,
              help = "plink MAP file", metavar = "MAP file"),
  make_option(c("-W", "--window"), default=15,
              help = "Sliding window size", metavar = "Window size"),
  make_option(c("-T", "--threshold"), default=0.05,
              help = "threshold of overlapping windows", metavar = "overlapping threshold"),
  make_option(c("-S", "--minSNP"), default=3,
              help = "minimum number of SNP", metavar = "minimum SNP"),
  make_option(c("-O", "--maxhom"), default=1,
              help = "maximum number of homozygous SNPs allowed", metavar = "max hom"),
  make_option(c("-L", "--maxmissing"), default=1,
              help = "maximum number of missing SNPs allowed", metavar = "max miss"),
  make_option(c("-G", "--maxgap"), default=1000000,
              help = "maximum gap between SNPs", metavar = "max gap"),
  make_option(c("-B", "--minlength"), default=1000,
              help = "minimum length of a window (BP)", metavar = "minimum length"),
  make_option(c("-D", "--mindensity"), default=0.1,
              help = "minimum density of SNPs (SNP/KB)", metavar = "minimum density"),
  make_option(c("-F", "--maxhomrun"), default=NULL,
              help = "maximum number of homozygous SNP allowed in whole run", metavar = "max hom in run"),
  make_option(c("-R", "--maxmissingrun"), default=NULL,
              help = "maximum number of missing SNP allowed in whole run", metavar = "max miss in run"),
  make_option(c("-C", "--chromosome"), default=NULL,
              help = "name of chromosome", metavar = "chromosome name")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

if(is.null(opt$ped) | is.null(opt$map) | is.null(opt$chromosome)) {
  print_help(opt_parser)
  stop("Please provide the ped and map files and chromosome name ", call. = FALSE)

}

genotypeFilePath <- opt$ped
mapFilePath <- opt$map
chrom <- gsub("[^[:alnum:]]", "", opt$chromosome)

slidingRUNS <- slidingRUNS.run(genotypeFile=opt$ped,
                               mapFile=opt$map,
                               windowSize=opt$window,
                               threshold=opt$threshold,
                               minSNP=opt$minSNP,
                               ROHet=TRUE,
                               maxOppWindow=opt$maxhom,
                               maxMissWindow=opt$maxmissing,
                               maxGap=opt$maxgap,
                               minLengthBps=opt$minlength,
                               minDensity=opt$mindensity,
                               maxOppRun=opt$maxhomrun,
                               maxMissRun=opt$maxmissingrun
                               )

fix_df <- slidingRUNS %>%
  rename(
    FID = group,
    IID = id,
    CHR = chrom,
    POS1 = from,
    POS2 = to,
    KB = lengthBps,
    NSNP = nSNP
    ) %>%
  mutate(
    PHE = -9.000,
    SNP1 = ".",
    SNP2 = ".",
    KB = KB / 1000,
    DENSITY = KB / NSNP,
    PHOM = ".",
    PHET = "."
    ) %>%
  relocate(NSNP, .after = KB) %>%
  relocate(PHE, .after = IID) %>%
  relocate(SNP1, .after = CHR) %>%
  relocate(SNP2, .after = SNP1)


write.table(fix_df, file = paste0(chrom, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

