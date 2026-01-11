#!/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(zoo))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

option_list = list(
  make_option(c("-T", "--tsv"), default=NULL,
              help = "list of tsv files", metavar = "tsv files"),
  make_option(c("-V", "--version"), default=NULL,
              help = "Nextflow process version", metavar = "process version")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser)

if(is.null(opt$tsv)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call. = FALSE)

}

ROHfile <- read.table(opt$tsv, sep = "\t", fill = TRUE, header = TRUE) %>%
  arrange(as.integer(sub("ROH_", "", ROH_ID)))

write.table(ROHfile, file = "combined.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

num_genes <- ROHfile %>%
  mutate(gene_count = str_count(genes, " ")) %>%   # count spaces
  select(ROH_ID, gene_count) %>%
  mutate(gene_count = replace_na(gene_count, 0)) %>%  # optional
  mutate(gene_count = as.integer(gene_count)) %>%
  {
    full_range <- min(.$gene_count):max(.$gene_count)

    data.frame(
      gene_count = full_range,
      n = as.vector(table(factor(.$gene_count, levels = full_range)))
    )
  }

ggplot(num_genes, aes(x=as.factor(gene_count), y=as.integer(n), fill = as.factor(gene_count))) +
  geom_bar(stat="identity") +
  xlab("number of Genes") +
  ylab("number of ROH") +
  scale_y_sqrt(
    limits = c(0, max(num_genes$n)),
    expand = c(0, 0)
  ) +
  theme(legend.position = "none") +
  theme_minimal()

ggsave(
  "genes_roh_mqc.png",
  width = 18,   # wide width in inches
  height = 6,   # reasonable height
  dpi = 300     # good resolution for print/publication
)

ROH_info <- ROHfile %>%
  select(ROH_ID, KB, NSNP, DENSITY)

k <- max(1, round(0.01 * nrow(ROH_info)))
ROH_info$y1_kb <- rollmean(ROH_info$KB, k = k, fill = NA, align = "center")
ROH_info$y1_nsnp <- rollmean(ROH_info$NSNP, k = k, fill = NA, align = "center")
ROH_info$y1_dens <- rollmean(1/ROH_info$DENSITY, k = k, fill = NA, align = "center")
ROH_info$ROH_num <- as.numeric(sub("ROH_", "", ROH_info$ROH_ID))
x_labels <- x_labels <- c(ROH_info$ROH_num[1], ROH_info$ROH_num[nrow(ROH_info)])

p1 <- ggplot(ROH_info, aes(x = ROH_num)) +
  geom_line(aes(y = y1_kb, group = 1), color = "blue") +
  scale_y_sqrt() +  # keeps zero values visible
  scale_x_continuous(
    breaks = x_labels,
    labels = c(ROH_info$ROH_ID[1], ROH_info$ROH_ID[nrow(ROH_info)])
  ) +
  labs(y = paste0("length ROH (KB) (grouped k = ",  k, ")"), x = "ROH ID") +
  theme_minimal()

p2 <- ggplot(ROH_info, aes(x = ROH_num)) +
  geom_line(aes(y = y1_nsnp, group = 1), color = "red") +
  scale_y_sqrt() +  # keeps zero values visible
  scale_x_continuous(
    breaks = x_labels,
    labels = c(ROH_info$ROH_ID[1], ROH_info$ROH_ID[nrow(ROH_info)])
  ) +
  labs(y = paste0("NSNP (grouped k = ",  k, ")"), x = "ROH ID") +
  theme_minimal()

p3 <- ggplot(ROH_info, aes(x = ROH_num)) +
  geom_line(aes(y = y1_dens, group = 1), color = "#187118") +
  scale_y_sqrt() +  # keeps zero values visible
  scale_x_continuous(
    breaks = x_labels,
    labels = c(ROH_info$ROH_ID[1], ROH_info$ROH_ID[nrow(ROH_info)])
  ) +
  labs(y = paste0("SNP density (SNP/KB) (grouped k = ",  k, ")"), x = "ROH ID") +
  theme_minimal()

combined <- p1 + p2 + p3
ggsave(
  "three_plots_mqc.png",
  combined,
  width = 18,   # wide width in inches
  height = 6,   # reasonable height
  dpi = 300     # good resolution for print/publication
)

num_samples <- ROHfile %>%
  mutate(sample_count = str_count(samples, ",") + 1) %>%   # count spaces
  select(ROH_ID, sample_count, samples) %>%
  arrange(-sample_count)

write.table(num_samples, file = "ROH_islands.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

count_samples <- num_samples %>%
  count(sample_count)

ggplot(count_samples, aes(x = sample_count, y = n, fill = as.factor(sample_count))) +
  geom_bar(stat = "identity") +
  scale_y_sqrt() +
  xlab("sample count") +
  ylab("number of ROH") +
  theme(legend.position = "none")

ggsave(
  "ROH_sample_mqc.png",
  width = 12,   # wide width in inches
  height = 6,   # reasonable height
  dpi = 300     # good resolution for print/publication
)

ver_zoo <- as.character(packageVersion("zoo"))
ver_patchwork <- as.character(packageVersion("patchwork"))
ver_R <- strsplit(R.version.string, " ")[[1]][3]

settings_list <- list(
  `Process version` = list(
    plots = opt$version
  ),
  `Tool version` = list(
    zoo =  ver_zoo,
    R = ver_R,
    patchwork = ver_patchwork
  )
)

write_yaml(settings_list, file = "versions.yml")