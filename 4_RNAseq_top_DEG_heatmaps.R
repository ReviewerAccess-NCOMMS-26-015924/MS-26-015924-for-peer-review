# =============================================================================
# 4_RNAseq_top_DEG_heatmaps.R
#
# Description:
# This script identifies the top differentially expressed genes independently
# at 1 and 4 weeks using limma, visualizes their expression in heatmaps, and
# exports genes from each row cluster for enrichment analysis.
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
  library(pheatmap)
  library(gridExtra)
})

# =============================================================================
# File paths and analysis parameters
# =============================================================================

normalized_file <- "astrocyte_normalized.RData"

top_n <- 500
n_row_clusters <- 4

if (!file.exists(normalized_file)) {
  stop("Input file not found: ", normalized_file, call. = FALSE)
}

load(normalized_file)

# Define the output path after load() so that a save_dir object stored in an
# older RData file cannot overwrite the portable path used by this script.
save_dir <- "output"

required_objects <- c("logtpm_mat", "col_meta", "row_meta")
missing_objects <- setdiff(required_objects, ls())

if (length(missing_objects) > 0) {
  stop(
    "Missing objects in the input file: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

group_order_1wk <- c(
  "astro_1wk_ETP0",
  "astro_1wk_ETP10",
  "astro_1wk_ETP10_CB"
)

group_order_4wk <- c(
  "astro_4wk_ETP0",
  "astro_4wk_ETP10",
  "astro_4wk_ETP10_CB"
)

group_all <- c(group_order_1wk, group_order_4wk)

# =============================================================================
# Validate and align sample metadata
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

missing_matrix_samples <- setdiff(col_meta$sample, colnames(logtpm_mat))
if (length(missing_matrix_samples) > 0) {
  stop(
    "Samples in col_meta are missing from logtpm_mat: ",
    paste(missing_matrix_samples, collapse = ", "),
    call. = FALSE
  )
}

col_meta_all <- col_meta[col_meta$group %in% group_all, , drop = FALSE]

missing_groups <- setdiff(group_all, unique(col_meta_all$group))
if (length(missing_groups) > 0) {
  stop(
    "Missing experimental groups: ",
    paste(missing_groups, collapse = ", "),
    call. = FALSE
  )
}

mat_all_samples <- logtpm_mat[
  ,
  col_meta_all$sample,
  drop = FALSE
]

# =============================================================================
# Select top genes independently at each time point
# =============================================================================

select_top_genes <- function(group_order, number_of_genes) {
  meta_subset <- col_meta_all[
    col_meta_all$group %in% group_order,
    ,
    drop = FALSE
  ]
  
  meta_subset$group <- factor(meta_subset$group, levels = group_order)
  expression_matrix <- logtpm_mat[
    ,
    meta_subset$sample,
    drop = FALSE
  ]
  
  design <- model.matrix(~ 0 + group, data = meta_subset)
  colnames(design) <- levels(meta_subset$group)
  
  fit <- lmFit(expression_matrix, design)
  
  contrast_matrix <- makeContrasts(
    contrasts = c(
      paste(group_order[1], "-", group_order[2]),
      paste(group_order[1], "-", group_order[3]),
      paste(group_order[2], "-", group_order[3])
    ),
    levels = design
  )
  
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  
  results <- topTable(fit, number = Inf, sort.by = "F")
  head(rownames(results), n = min(number_of_genes, nrow(results)))
}

top_ids_1wk <- select_top_genes(group_order_1wk, top_n)
top_ids_4wk <- select_top_genes(group_order_4wk, top_n)

# =============================================================================
# Calculate gene-wise Z-scores using all 1- and 4-week samples
# =============================================================================

scale_using_all_samples <- function(gene_ids, target_samples) {
  expression_all <- mat_all_samples[gene_ids, , drop = FALSE]
  scaled_all <- t(scale(t(expression_all)))
  
  # Genes with zero variance across all samples cannot be Z-score transformed.
  keep_finite <- apply(scaled_all, 1, function(x) all(is.finite(x)))
  scaled_target <- scaled_all[keep_finite, target_samples, drop = FALSE]
  
  # Correlation-based row clustering requires variation within the time point.
  keep_variable <- apply(scaled_target, 1, sd) > 0
  scaled_target <- scaled_target[keep_variable, , drop = FALSE]
  
  if (nrow(scaled_target) < 2) {
    stop("Too few variable genes remained for heatmap clustering.", call. = FALSE)
  }
  
  scaled_target
}

samples_1wk <- col_meta_all$sample[
  col_meta_all$group %in% group_order_1wk
]
samples_4wk <- col_meta_all$sample[
  col_meta_all$group %in% group_order_4wk
]

mat_1wk <- scale_using_all_samples(top_ids_1wk, samples_1wk)
mat_4wk <- scale_using_all_samples(top_ids_4wk, samples_4wk)

# =============================================================================
# Order samples by group and by similarity within each group
# =============================================================================

order_within_groups <- function(scaled_matrix, group_order) {
  metadata <- col_meta_all[
    match(colnames(scaled_matrix), col_meta_all$sample),
    ,
    drop = FALSE
  ]
  metadata$group <- factor(metadata$group, levels = group_order)
  
  ordered_samples <- character(0)
  
  for (group_name in group_order) {
    group_samples <- metadata$sample[metadata$group == group_name]
    
    if (length(group_samples) > 2) {
      sample_distance <- dist(t(scaled_matrix[, group_samples, drop = FALSE]))
      sample_tree <- hclust(sample_distance)
      group_samples <- group_samples[sample_tree$order]
    }
    
    ordered_samples <- c(ordered_samples, group_samples)
  }
  
  ordered_samples
}

sample_order_1wk <- order_within_groups(mat_1wk, group_order_1wk)
sample_order_4wk <- order_within_groups(mat_4wk, group_order_4wk)

mat_1wk <- mat_1wk[, sample_order_1wk, drop = FALSE]
mat_4wk <- mat_4wk[, sample_order_4wk, drop = FALSE]

# =============================================================================
# Define a common heatmap color scale
# =============================================================================

global_max <- ceiling(
  max(abs(c(mat_1wk, mat_4wk)), na.rm = TRUE)
)

if (!is.finite(global_max) || global_max <= 0) {
  stop("Could not define the heatmap color scale.", call. = FALSE)
}

heatmap_breaks <- seq(-global_max, global_max, length.out = 101)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# =============================================================================
# Column annotations
# =============================================================================

make_annotation <- function(samples, group_order) {
  sample_groups <- col_meta_all$group[
    match(samples, col_meta_all$sample)
  ]
  
  data.frame(
    group = factor(sample_groups, levels = group_order),
    row.names = samples
  )
}

annotation_1wk <- make_annotation(colnames(mat_1wk), group_order_1wk)
annotation_4wk <- make_annotation(colnames(mat_4wk), group_order_4wk)

annotation_colors_1wk <- list(
  group = setNames(c("lightgray", "red", "pink"), group_order_1wk)
)

annotation_colors_4wk <- list(
  group = setNames(c("lightgray", "red", "pink"), group_order_4wk)
)

# =============================================================================
# Draw heatmaps
# =============================================================================

heatmap_1wk <- pheatmap(
  mat_1wk,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  annotation_col = annotation_1wk,
  annotation_colors = annotation_colors_1wk,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 7,
  treeheight_row = 50,
  main = "1wk",
  silent = TRUE
)

heatmap_4wk <- pheatmap(
  mat_4wk,
  cluster_rows = TRUE,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  annotation_col = annotation_4wk,
  annotation_colors = annotation_colors_4wk,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 7,
  treeheight_row = 50,
  main = "4wk",
  silent = TRUE
)

# =============================================================================
# Save the combined heatmap
# =============================================================================

heatmap_file <- file.path(
  save_dir,
  paste0("heatmap_top", top_n, "_genes_1wk_4wk_global_zscore.pdf")
)

pdf(heatmap_file, width = 10, height = 8)
grid.arrange(heatmap_1wk$gtable, heatmap_4wk$gtable, ncol = 2)
dev.off()

message("Saved: ", heatmap_file)

# =============================================================================
# Export genes from each row cluster for enrichment analysis
# =============================================================================

cluster_output_dir <- file.path(save_dir, "cluster_genes_for_Enrichr")
if (!dir.exists(cluster_output_dir)) {
  dir.create(cluster_output_dir, recursive = TRUE)
}

ensembl_to_symbol <- setNames(row_meta$SYMBOL, row_meta$ENSEMBL)

export_cluster_genes <- function(
    scaled_matrix,
    row_tree,
    timepoint_label,
    number_of_clusters) {
  
  if (nrow(scaled_matrix) < number_of_clusters) {
    stop(
      "The number of row clusters exceeds the number of displayed genes.",
      call. = FALSE
    )
  }
  
  cluster_assignment <- cutree(row_tree, k = number_of_clusters)
  
  for (cluster_id in sort(unique(cluster_assignment))) {
    ensembl_ids <- names(cluster_assignment)[
      cluster_assignment == cluster_id
    ]
    
    gene_symbols <- unname(ensembl_to_symbol[ensembl_ids])
    gene_symbols <- unique(gene_symbols[
      !is.na(gene_symbols) & nzchar(gene_symbols)
    ])
    
    output_file <- file.path(
      cluster_output_dir,
      paste0(
        timepoint_label,
        "_cluster_",
        cluster_id,
        "_genes.csv"
      )
    )
    
    write.csv(
      data.frame(gene = gene_symbols),
      file = output_file,
      row.names = FALSE,
      quote = FALSE
    )
    
    message(
      "Saved: ", output_file,
      " (", length(gene_symbols), " unique gene symbols)"
    )
  }
}

export_cluster_genes(
  mat_1wk,
  heatmap_1wk$tree_row,
  "1wk",
  n_row_clusters
)

export_cluster_genes(
  mat_4wk,
  heatmap_4wk$tree_row,
  "4wk",
  n_row_clusters
)
