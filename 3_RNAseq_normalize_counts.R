# =============================================================================
# 3_RNAseq_normalize_counts.R
#
# Description:
# This script imports the gene-level featureCounts matrix, assigns sample
# metadata, filters lowly expressed genes with edgeR, calculates TMM
# normalization factors, and generates TPM and log2(TPM + 1) matrices.
#
# Author:
# Momoka Hikosaka
# Kyoto University
#
# Associated publication:
# Hikosaka et al. in review
# "Glial metabolic states in human sporadic amyotrophic lateral sclerosis
# uncovered by spatial single-cell proteomics."
#
# Input:
# - gene_counts.txt
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# =============================================================================
# File paths
# =============================================================================
counts_file <- "gene_counts.txt"
save_dir <- "output"

if (!file.exists(counts_file)) {
  stop("Count file not found: ", counts_file, call. = FALSE)
}

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# =============================================================================
# Sample metadata
# =============================================================================

sample_meta <- data.frame(
  sample = c(
    "astro1_1", "astro1_2", "astro1_3",
    "astro1_4", "astro1_5", "astro1_6",
    "astro1_7", "astro1_8", "astro1_9",
    "astro4_1", "astro4_2", "astro4_3",
    "astro4_4", "astro4_5", "astro4_6",
    "astro4_7", "astro4_8", "astro4_9"
  ),
  timepoint = rep(c("1wk", "4wk"), each = 9),
  ETP_treat = rep(
    c(rep("ETP0", 3), rep("ETP10", 6)),
    times = 2
  ),
  reagent = rep(
    c(rep("-", 6), rep("CB", 3)),
    times = 2
  ),
  replicate = rep(1:3, times = 6),
  stringsAsFactors = FALSE
)

# =============================================================================
# Import the featureCounts matrix
# =============================================================================

base_data <- read.delim(
  counts_file,
  header = TRUE,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length"
)
missing_columns <- setdiff(required_columns, colnames(base_data))

if (length(missing_columns) > 0) {
  stop(
    "Missing featureCounts columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(base_data$Geneid)) {
  stop("Duplicated Geneid values were detected.", call. = FALSE)
}

gene_length_bp <- setNames(base_data$Length, base_data$Geneid)

count_columns <- setdiff(colnames(base_data), required_columns)
if (length(count_columns) == 0) {
  stop("No sample count columns were found.", call. = FALSE)
}

counts_matrix <- as.matrix(base_data[, count_columns, drop = FALSE])
storage.mode(counts_matrix) <- "numeric"
rownames(counts_matrix) <- base_data$Geneid

if (anyNA(counts_matrix) || any(counts_matrix < 0)) {
  stop("The count matrix contains missing or negative values.", call. = FALSE)
}

# featureCounts uses BAM paths as column names. Retain only the sample name.
sample_list <- basename(gsub("\\\\", "/", colnames(counts_matrix)))
sample_list <- sub("\\.sorted\\.bam$", "", sample_list)

if (anyDuplicated(sample_list)) {
  stop("Sample names derived from BAM filenames are not unique.", call. = FALSE)
}

unmatched_samples <- setdiff(sample_list, sample_meta$sample)
missing_samples <- setdiff(sample_meta$sample, sample_list)

if (length(unmatched_samples) > 0 || length(missing_samples) > 0) {
  stop(
    paste0(
      "Sample names do not match the metadata. ",
      "Unmatched count columns: ",
      paste(unmatched_samples, collapse = ", "),
      "; missing expected samples: ",
      paste(missing_samples, collapse = ", ")
    ),
    call. = FALSE
  )
}

colnames(counts_matrix) <- sample_list
sample_meta <- sample_meta[match(sample_list, sample_meta$sample), , drop = FALSE]

if (anyNA(sample_meta$sample)) {
  stop("Sample metadata could not be aligned to the count matrix.", call. = FALSE)
}

# =============================================================================
# Define experimental groups
# =============================================================================

timepoint <- factor(sample_meta$timepoint, levels = c("1wk", "4wk"))
ETP_treat <- factor(sample_meta$ETP_treat, levels = c("ETP0", "ETP10"))
reagent <- factor(sample_meta$reagent, levels = c("-", "CB"))
replicate <- sample_meta$replicate

treatment <- ifelse(
  reagent == "CB",
  paste0(ETP_treat, "_CB"),
  as.character(ETP_treat)
)

group_levels <- paste0(
  "astro_",
  c(
    "1wk_ETP0", "1wk_ETP10", "1wk_ETP10_CB",
    "4wk_ETP0", "4wk_ETP10", "4wk_ETP10_CB"
  )
)

group <- factor(
  paste("astro", timepoint, treatment, sep = "_"),
  levels = group_levels
)

if (anyNA(group)) {
  stop("One or more samples could not be assigned to a group.", call. = FALSE)
}

# =============================================================================
# edgeR filtering and TMM normalization
# =============================================================================

dge <- DGEList(counts = counts_matrix, group = group)

keep <- filterByExpr(
  dge,
  group = group,
  min.count = 1,
  min.total.count = 1
)

if (!any(keep)) {
  stop("No genes remained after expression filtering.", call. = FALSE)
}

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# =============================================================================
# TPM calculation
# =============================================================================

filtered_gene_length_bp <- gene_length_bp[rownames(dge$counts)]

if (anyNA(filtered_gene_length_bp) || any(filtered_gene_length_bp <= 0)) {
  stop("Missing or non-positive gene lengths were detected.", call. = FALSE)
}

gene_length_kb <- filtered_gene_length_bp / 1000
reads_per_kb <- sweep(dge$counts, 1, gene_length_kb, FUN = "/")
per_sample_scaling <- colSums(reads_per_kb) / 1e6

if (any(per_sample_scaling <= 0)) {
  stop("TPM scaling failed because a sample has no count signal.", call. = FALSE)
}

tpm_mat <- sweep(reads_per_kb, 2, per_sample_scaling, FUN = "/")
logtpm_mat <- log2(tpm_mat + 1)

# =============================================================================
# Gene annotation
# =============================================================================

ens_ids <- rownames(dge$counts)
ens_stripped <- sub("\\.\\d+$", "", ens_ids)

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = ens_stripped,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# =============================================================================
# Assemble output objects
# =============================================================================

col_meta <- data.frame(
  sample = colnames(dge$counts),
  timepoint = as.character(timepoint),
  ETP_treat = as.character(ETP_treat),
  reagent = as.character(reagent),
  replicate = replicate,
  group = as.character(group),
  lib.size = dge$samples$lib.size,
  norm.factor = dge$samples$norm.factors,
  stringsAsFactors = FALSE,
  row.names = colnames(dge$counts)
)

row_meta <- data.frame(
  ENSEMBL = ens_ids,
  ENSEMBL_stripped = ens_stripped,
  SYMBOL = unname(symbol_map[ens_stripped]),
  Length_bp = unname(filtered_gene_length_bp),
  stringsAsFactors = FALSE,
  row.names = ens_ids
)

expr_data <- list(
  assays = list(
    raw = dge$counts,
    tpm = tpm_mat,
    logTPM = logtpm_mat
  ),
  colData = col_meta,
  rowData = row_meta,
  group = group,
  norm.factors = dge$samples$norm.factors
)

# =============================================================================
# Save R objects and text tables
# =============================================================================

workspace_file <- file.path(save_dir, "astrocyte_normalized.RData")
tpm_file <- file.path(save_dir, "astrocyte_normalized_TPM.csv")
logtpm_file <- file.path(save_dir, "astrocyte_normalized_logTPM.csv")

save(
  dge,
  group,
  timepoint,
  ETP_treat,
  reagent,
  replicate,
  sample_list,
  sample_meta,
  col_meta,
  row_meta,
  tpm_mat,
  logtpm_mat,
  expr_data,
  file = workspace_file
)

tpm_out <- data.frame(
  Geneid = row_meta$ENSEMBL,
  Symbol = row_meta$SYMBOL,
  tpm_mat,
  check.names = FALSE,
  row.names = NULL
)

logtpm_out <- data.frame(
  Geneid = row_meta$ENSEMBL,
  Symbol = row_meta$SYMBOL,
  logtpm_mat,
  check.names = FALSE,
  row.names = NULL
)

write.csv(tpm_out, tpm_file, row.names = FALSE, na = "")
write.csv(logtpm_out, logtpm_file, row.names = FALSE, na = "")

message("Saved: ", workspace_file)
message("Saved: ", tpm_file)
message("Saved: ", logtpm_file)