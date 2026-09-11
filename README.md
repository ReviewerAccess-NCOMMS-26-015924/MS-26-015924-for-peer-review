## 4. iPSC-Derived Astrocyte Bulk RNA-seq Analysis (WSL and R)

This module contains the scripts used to preprocess and analyze bulk RNA-seq data obtained from iPSC-derived astrocytes at 1 and 4 weeks.

The workflow includes:

* Read preprocessing
* Alignment to the human reference genome
* Gene-level quantification
* Expression normalization
* Differential-expression analysis
* Heatmap visualization
* Preparation of pre-ranked gene lists for gene set enrichment analysis (GSEA)

The raw RNA-seq data have been deposited in the Gene Expression Omnibus (GEO) under accession number **GSE###**.

### `1_RNAseq_preprocessing.sh`

This script performs read preprocessing using fastp and paired-end alignment to the human reference genome using HISAT2. The resulting alignments are coordinate-sorted, indexed, and summarized using samtools.

### `2_RNAseq_featurecounts.sh`

This script quantifies gene-level expression from coordinate-sorted BAM files using featureCounts. All available samples are quantified in a single run to generate a combined count matrix and its corresponding read-assignment summary.

### `3_RNAseq_normalize_counts.R`

This script imports the gene-level featureCounts matrix, assigns sample metadata, filters lowly expressed genes using edgeR, calculates trimmed mean of M-values (TMM) normalization factors, and generates TPM and `log2(TPM + 1)` expression matrices.

### `4_RNAseq_top_DEG_heatmaps.R`

This script independently identifies the top differentially expressed genes at 1 and 4 weeks using limma, visualizes their expression in heatmaps, and exports the genes assigned to each heatmap row cluster for enrichment analysis.

### `5_RNAseq_prepare_GSEA_ranked_list.R`

This script uses limma to calculate moderated *t*-statistics for the comparison between ETP10- and ETP0-treated 1-week astrocytes and exports a pre-ranked gene list for GSEA Pre-ranked analysis.
