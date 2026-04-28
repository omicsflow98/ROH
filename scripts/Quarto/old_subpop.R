#!/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(plotly))
suppressMessages(library(DT))


nroh_subpop <- function(subpopulation) {
  if (!is.null(subpopulation)) {
    return(nrow(subpopulation)) 
  } else {
    return (0)
  }
}

subpopulations <- read.table(subpop, 
                             sep = "\t", 
                             fill = TRUE, 
                             header = TRUE,
                             check.names = FALSE)
colnames(subpopulations) <- gsub("\\s+", "_", colnames(subpopulations))

if (hom) {
  ROHfile <- read.table(homfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    arrange(as.integer(sub("ROH_", "", ROH_ID)))
  
  ROH_subpop <- ROHfile %>%
    select(ROH_ID, samples) %>%
    separate_wider_delim(
      samples,
      delim = ",",
      names_sep = "_",
      too_few = "align_start"
    )
  
  # Convert subpopulations to matrix
  mat <- as.matrix(subpopulations[ , -1])
  rownames(mat) <- subpopulations$SampleID
  
  # Identify ID columns in ROH_subpop
  id_cols <- grep("^samples_", colnames(ROH_subpop), value = TRUE)
  
  # Shared key column name
  key_col <- "ROH_ID"   # <-- change to your actual column name
  
  # Initialize empty list
  column_groups <- setNames(vector("list", ncol(mat)), colnames(mat))
  
  # Loop over ROH_subpop rows
  for (i in seq_len(nrow(ROH_subpop))) {
    
    ids <- unlist(ROH_subpop[i, id_cols])
    ids <- ids[!is.na(ids)]
    ids <- intersect(ids, rownames(mat))
    if (length(ids) == 0) next
    
    sub_mat <- mat[ids, , drop = FALSE]
    matching_cols <- colnames(sub_mat)[colSums(sub_mat) > 0]
    
    # Get shared ID from ROH_subpop
    shared_value <- ROH_subpop[[key_col]][i]
    
    for (col in matching_cols) {
      column_groups[[col]] <- c(column_groups[[col]], shared_value)
    }
  }
  
  # Now subset ROHfile instead of ROH_subpop
  ROHfile_subsets <- lapply(column_groups, function(keys) {
    if (length(keys) > 0) {
      ROHfile %>% filter(.data[[key_col]] %in% unique(keys))
    } else {
      NULL
    }
  })
  
  rm(mat, column_groups, sub_mat, col, i, id_cols, ids, key_col, matching_cols, shared_value)
  
  #number of samples per subpopulation
  num_samples <- colSums(subpopulations[,-1])
  num_samples_df <- enframe(num_samples) %>%
    rename(Subpopulation = name, `Number of Samples` = value)
  num_samples_df$Subpopulation <- gsub("_", " ", num_samples_df$Subpopulation)
  
  #number of ROH in each subpopulation
  rows_subpop_ROH <- sapply(ROHfile_subsets, nroh_subpop)
  rows_df_ROH <- enframe(rows_subpop_ROH)
  rows_df_ROH$name <- gsub("_", " ", rows_df_ROH$name)
  
  #ROH Boxplot
  for (nm in names(ROHfile_subsets)) {
    ROHfile_subsets[[nm]]$subpopulation <- nm
  }
  joint_ROHfile <- bind_rows(ROHfile_subsets) %>%
    select(subpopulation, KB) %>%
    mutate(KB = replace_na(KB, 0))
  
  joint_ROHfile$subpopulation <- gsub("_", " ", joint_ROHfile$subpopulation)
}

if (het) {
  HRRfile <- read.table(hetfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    arrange(as.integer(sub("HRR_", "", HRR_ID)))
  
  HRR_subpop <- HRRfile %>%
    select(HRR_ID, samples) %>%
    separate_wider_delim(
      samples,
      delim = ",",
      names_sep = "_",
      too_few = "align_start"
    )
  
  # Convert subpopulations to matrix
  mat <- as.matrix(subpopulations[ , -1])
  rownames(mat) <- subpopulations$SampleID
  
  # Identify ID columns in ROH_subpop
  id_cols <- grep("^samples_", colnames(HRR_subpop), value = TRUE)
  
  # Shared key column name
  key_col <- "HRR_ID"   # <-- change to your actual column name
  
  # Initialize empty list
  column_groups <- setNames(vector("list", ncol(mat)), colnames(mat))
  
  # Loop over ROH_subpop rows
  for (i in seq_len(nrow(HRR_subpop))) {
    
    ids <- unlist(HRR_subpop[i, id_cols])
    ids <- ids[!is.na(ids)]
    ids <- intersect(ids, rownames(mat))
    if (length(ids) == 0) next
    
    sub_mat <- mat[ids, , drop = FALSE]
    matching_cols <- colnames(sub_mat)[colSums(sub_mat) > 0]
    
    # Get shared ID from ROH_subpop
    shared_value <- HRR_subpop[[key_col]][i]
    
    for (col in matching_cols) {
      column_groups[[col]] <- c(column_groups[[col]], shared_value)
    }
  }
  
  # Now subset ROHfile instead of ROH_subpop
  HRRfile_subsets <- lapply(column_groups, function(keys) {
    if (length(keys) > 0) {
      HRRfile %>% filter(.data[[key_col]] %in% unique(keys))
    } else {
      NULL
    }
  })
  
  rm(mat, column_groups, sub_mat, col, i, id_cols, ids, key_col, matching_cols, shared_value)
  
  #number of HRR in each subpopulation
  rows_subpop_HRR <- sapply(HRRfile_subsets, nroh_subpop)
  rows_df_HRR <- enframe(rows_subpop_HRR)
  rows_df_HRR$name <- gsub("_", " ", rows_df_HRR$name)
  
  #ROH Boxplot
  for (nm in names(HRRfile_subsets)) {
    HRRfile_subsets[[nm]]$subpopulation <- nm
  }
  joint_HRRfile <- bind_rows(HRRfile_subsets) %>%
    select(subpopulation, KB) %>%
    mutate(KB = replace_na(KB, 0))
  
  joint_HRRfile$subpopulation <- gsub("_", " ", joint_HRRfile$subpopulation)
}
