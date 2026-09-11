#!/bin/bash
# =============================================================================
# 2_RNAseq_featurecounts.sh
#
# Description:
# This script quantifies gene-level expression from coordinate-sorted BAM files
# using featureCounts. All available samples are quantified in a single run to
# produce one count matrix and its assignment summary.
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
# Requirements:
# - Bash 4 or later
# - featureCounts (Subread)
# - Coordinate-sorted BAM files 
# - GTF annotation 
#
# An optional third argument specifies the output directory:
#   bash 03_run_featureCounts.sh BAM_DIR GTF_FILE COUNTS_DIR
#
# The paths may alternatively be supplied through the BAM_DIR, GTF_FILE, and
# COUNTS_DIR environment variables.
# =============================================================================

set -euo pipefail

# =============================================================================
# Paths
# =============================================================================

BAM_DIR="${1:-${BAM_DIR:-}}"
GTF_FILE="${2:-${GTF_FILE:-}}"
COUNTS_DIR="${3:-${COUNTS_DIR:-}}"

if [[ -z "$BAM_DIR" || -z "$GTF_FILE" ]]; then
  echo "Usage: bash $0 BAM_DIR GTF_FILE [COUNTS_DIR]" >&2
  exit 2
fi

if [[ -z "$COUNTS_DIR" ]]; then
  COUNTS_DIR="${BAM_DIR}/featureCounts"
fi

COUNTS_OUT="${COUNTS_DIR}/gene_counts.txt"
VERSION_FILE="${COUNTS_DIR}/software_versions.txt"

# =============================================================================
# Analysis parameters
# =============================================================================

THREADS=16

# Library strandedness:
#   0 = unstranded
#   1 = stranded
#   2 = reversely stranded
STRANDEDNESS=0

# =============================================================================
# Helper functions
# =============================================================================

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "ERROR: Required command not found: $command_name" >&2
    echo "Install Subread in the analysis environment before running this script." >&2
    exit 127
  fi
}

report_error() {
  local exit_code=$?
  local line_number="$1"
  echo "ERROR: Pipeline terminated at line ${line_number} (exit code ${exit_code})." >&2
  exit "$exit_code"
}

trap 'report_error "$LINENO"' ERR

# =============================================================================
# Preflight checks
# =============================================================================

require_command featureCounts

if [[ ! -d "$BAM_DIR" ]]; then
  echo "ERROR: BAM directory not found: $BAM_DIR" >&2
  exit 1
fi

if [[ ! -s "$GTF_FILE" ]]; then
  echo "ERROR: GTF annotation not found or empty: $GTF_FILE" >&2
  exit 1
fi

if [[ "$STRANDEDNESS" != "0" && "$STRANDEDNESS" != "1" && "$STRANDEDNESS" != "2" ]]; then
  echo "ERROR: STRANDEDNESS must be 0, 1, or 2." >&2
  exit 1
fi

mkdir -p "$COUNTS_DIR"

bam_files=("$BAM_DIR"/*.sorted.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
  echo "ERROR: No '*.sorted.bam' files found in: $BAM_DIR" >&2
  exit 1
fi

for bam_file in "${bam_files[@]}"; do
  if [[ ! -s "$bam_file" ]]; then
    echo "ERROR: BAM file is empty: $bam_file" >&2
    exit 1
  fi
done

{
  echo "Analysis started: $(date --iso-8601=seconds)"
  featureCounts -v 2>&1
} > "$VERSION_FILE"

echo "Input BAM files (${#bam_files[@]}):"
printf '  %s\n' "${bam_files[@]}"
echo

# =============================================================================
# Gene-level quantification
# =============================================================================

echo "Running featureCounts"

featureCounts \
  -T "$THREADS" \
  -p \
  --countReadPairs \
  -s "$STRANDEDNESS" \
  -a "$GTF_FILE" \
  -o "$COUNTS_OUT" \
  "${bam_files[@]}"

{
  echo "Analysis finished: $(date --iso-8601=seconds)"
  echo "Input BAM files: ${#bam_files[@]}"
  echo "Strandedness: $STRANDEDNESS"
  echo "Count matrix: $COUNTS_OUT"
  echo "Assignment summary: ${COUNTS_OUT}.summary"
} > "${COUNTS_DIR}/featureCounts_run_summary.txt"

echo
echo "featureCounts completed successfully."
echo "Count matrix: $COUNTS_OUT"
echo "Assignment summary: ${COUNTS_OUT}.summary"