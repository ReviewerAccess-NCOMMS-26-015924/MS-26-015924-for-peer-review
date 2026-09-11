# =============================================================================
# 5_RNAseq_prepare_GSEA_ranked_list.R
#
# Description:
# This script uses limma to calculate moderated t-statistics for the comparison
# between ETP10 and ETP0 in 1-week astrocytes and exports a pre-ranked gene list 
# for GSEA.
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
# - astrocyte_normalized.RData
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(limma)
})

# =============================================================================
# File paths
# =============================================================================

normalized_file <- "astrocyte_normalized.RData"

if (!file.exists(normalized_file)) {
  stop("Input file not found: ", normalized_file, call. = FALSE)
}

load(normalized_file)

output_dir <- file.path("output", "GSEA_input")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

required_objects <- c("logtpm_mat", "col_meta", "row_meta")
missing_objects <- setdiff(required_objects, ls())

if (length(missing_objects) > 0) {
  stop(
    "Missing objects in the input file: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# =============================================================================
# Define the comparison
# =============================================================================

group_order <- c(
  "astro_1wk_ETP10",
  "astro_1wk_ETP0"
)

contrast_name <- "astro_1wk_ETP10_vs_ETP0"

# The contrast is ETP10 minus ETP0. Therefore, positive ranking
# scores indicate higher expression in ETP10, whereas negative scores indicate
# higher expression in ETP0.

# =============================================================================
# Validate and select samples
# =============================================================================

required_metadata <- c("sample", "group")
missing_metadata <- setdiff(required_metadata, colnames(col_meta))

if (length(missing_metadata) > 0) {
  stop(
    "Missing col_meta columns: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(col_meta$sample)) {
  stop("Duplicated sample names were detected in col_meta.", call. = FALSE)
}

metadata_subset <- col_meta[
  col_meta$group %in% group_order,
  ,
  drop = FALSE
]

missing_groups <- setdiff(group_order, unique(metadata_subset$group))
if (length(missing_groups) > 0) {
  stop(
    "Missing experimental groups: ",
    paste(missing_groups, collapse = ", "),
    call. = FALSE
  )
}

group_counts <- table(metadata_subset$group)
if (any(group_counts[group_order] < 2)) {
  stop("At least two samples per group are required.", call. = FALSE)
}

metadata_subset$group <- factor(
  metadata_subset$group,
  levels = group_order
)
metadata_subset <- metadata_subset[
  order(metadata_subset$group, metadata_subset$sample),
  ,
  drop = FALSE
]

missing_samples <- setdiff(metadata_subset$sample, colnames(logtpm_mat))
if (length(missing_samples) > 0) {
  stop(
    "Samples missing from logtpm_mat: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}

expression_matrix <- logtpm_mat[
  ,
  metadata_subset$sample,
  drop = FALSE
]

# =============================================================================
# Calculate moderated t-statistics with limma
# =============================================================================

design <- model.matrix(~ 0 + group, data = metadata_subset)
colnames(design) <- levels(metadata_subset$group)

fit <- lmFit(expression_matrix, design)

contrast_matrix <- makeContrasts(
  astro_1wk_ETP10 - astro_1wk_ETP0,
  levels = design
)
colnames(contrast_matrix) <- contrast_name

fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

results <- topTable(
  fit,
  coef = contrast_name,
  number = Inf,
  sort.by = "none"
)
results$ENSEMBL <- rownames(results)

# =============================================================================
# Add gene symbols and resolve duplicated symbols
# =============================================================================

required_row_metadata <- c("ENSEMBL", "SYMBOL")
missing_row_metadata <- setdiff(required_row_metadata, colnames(row_meta))

if (length(missing_row_metadata) > 0) {
  stop(
    "Missing row_meta columns: ",
    paste(missing_row_metadata, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(row_meta$ENSEMBL)) {
  stop("Duplicated ENSEMBL identifiers were detected in row_meta.", call. = FALSE)
}

symbol_lookup <- setNames(row_meta$SYMBOL, row_meta$ENSEMBL)
results$SYMBOL <- unname(symbol_lookup[results$ENSEMBL])

results <- results[
  !is.na(results$SYMBOL) &
    nzchar(results$SYMBOL) &
    is.finite(results$t),
  ,
  drop = FALSE
]

# When multiple ENSEMBL identifiers map to the same symbol, retain the entry
# with the largest absolute moderated t-statistic.
results <- results[
  order(-abs(results$t), results$SYMBOL),
  ,
  drop = FALSE
]
results <- results[!duplicated(results$SYMBOL), , drop = FALSE]

# GSEA pre-ranked input is ordered from the largest to the smallest score.
results <- results[
  order(-results$t, results$SYMBOL),
  ,
  drop = FALSE
]

if (nrow(results) == 0) {
  stop("No genes remained for GSEA ranking.", call. = FALSE)
}

# =============================================================================
# Export the ranked list and complete statistics
# =============================================================================

ranked_list <- results[, c("SYMBOL", "t"), drop = FALSE]

ranked_file <- file.path(
  output_dir,
  paste0(contrast_name, ".rnk")
)

statistics_file <- file.path(
  output_dir,
  paste0(contrast_name, "_limma_results.tsv")
)

write.table(
  ranked_list,
  file = ranked_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  results[, c(
    "ENSEMBL", "SYMBOL", "logFC", "AveExpr",
    "t", "P.Value", "adj.P.Val", "B"
  )],
  file = statistics_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


message("Saved: ", ranked_file)
message("Saved: ", statistics_file)
