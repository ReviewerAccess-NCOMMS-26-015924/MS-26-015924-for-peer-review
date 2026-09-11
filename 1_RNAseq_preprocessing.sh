#!/bin/bash
# =============================================================================
# 1_RNAseq_preprocessing.sh
#
# Description:
# This script performs read preprocessing with fastp and paired-end alignment
# to the human reference genome with HISAT2. Alignments are coordinate-sorted,
# indexed, and summarized with samtools.
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
# - fastp
# - HISAT2
# - samtools
# - Reference files
#
# The three paths may alternatively be supplied through the INPUT_BASE,
# OUTPUT_DIR, and REF_DIR environment variables.
# =============================================================================

set -euo pipefail

# =============================================================================
# Paths
# =============================================================================

INPUT_BASE="${1:-${INPUT_BASE:-}}"
OUTPUT_DIR="${2:-${OUTPUT_DIR:-}}"
REF_DIR="${3:-${REF_DIR:-}}"

if [[ -z "$INPUT_BASE" || -z "$OUTPUT_DIR" || -z "$REF_DIR" ]]; then
  echo "Usage: bash $0 INPUT_BASE OUTPUT_DIR REF_DIR" >&2
  exit 2
fi

HISAT2_INDEX="${REF_DIR}/GRCh38_ens_hisat2"
SPLICE_SITES="${REF_DIR}/GRCh38_ens.splicesites.txt"

# =============================================================================
# Computational parameters
# =============================================================================

FASTP_THREADS=4
HISAT2_THREADS=16
SAMTOOLS_THREADS=8
SAMTOOLS_SORT_MEMORY="1G"

QUALIFIED_QUALITY_PHRED=15
LENGTH_REQUIRED=30

# =============================================================================
# Samples
# =============================================================================

SAMPLES=(
  astro1_1
  astro1_2
  astro1_3
  astro1_4
  astro1_5
  astro1_6
  astro1_7
  astro1_8
  astro1_9
  astro4_1
  astro4_2
  astro4_3
  astro4_4
  astro4_5
  astro4_6
  astro4_7
  astro4_8
  astro4_9
)

# =============================================================================
# Helper functions
# =============================================================================

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "ERROR: Required command not found: $command_name" >&2
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

require_command fastp
require_command hisat2
require_command samtools

if [[ ! -d "$INPUT_BASE" ]]; then
  echo "ERROR: Input directory not found: $INPUT_BASE" >&2
  exit 1
fi

if [[ ! -d "$REF_DIR" ]]; then
  echo "ERROR: Reference directory not found: $REF_DIR" >&2
  exit 1
fi

if [[ ! -e "${HISAT2_INDEX}.1.ht2" && ! -e "${HISAT2_INDEX}.1.ht2l" ]]; then
  echo "ERROR: HISAT2 index not found: $HISAT2_INDEX" >&2
  echo "Run 01_prepare_reference.sh before this script." >&2
  exit 1
fi

if [[ ! -s "$SPLICE_SITES" ]]; then
  echo "ERROR: Splice-site file not found or empty: $SPLICE_SITES" >&2
  echo "Run 01_prepare_reference.sh before this script." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Record the software versions used for the analysis.
{
  echo "Analysis started: $(date --iso-8601=seconds)"
  echo "fastp: $(fastp --version 2>&1 | head -n 1)"
  echo "HISAT2: $(hisat2 --version 2>&1 | head -n 1)"
  echo "samtools: $(samtools --version 2>&1 | head -n 1)"
} > "${OUTPUT_DIR}/software_versions.txt"

# =============================================================================
# Process each sample: fastp -> HISAT2 -> sorted BAM
# =============================================================================

processed_samples=0
skipped_samples=0

for sample in "${SAMPLES[@]}"; do
  echo "Processing sample: $sample"
  date

  sample_dir="${INPUT_BASE}/${sample}"

  if [[ ! -d "$sample_dir" ]]; then
    echo "WARNING: Sample directory not found; skipping: $sample_dir" >&2
    ((skipped_samples += 1))
    continue
  fi

  r1_candidates=("$sample_dir"/"${sample}"*_1.fq.gz)
  r2_candidates=("$sample_dir"/"${sample}"*_2.fq.gz)

  if [[ ${#r1_candidates[@]} -ne 1 || ${#r2_candidates[@]} -ne 1 ]]; then
    echo "WARNING: Expected exactly one R1 and one R2 file; skipping: $sample" >&2
    echo "  R1 candidates: ${#r1_candidates[@]}" >&2
    echo "  R2 candidates: ${#r2_candidates[@]}" >&2
    ((skipped_samples += 1))
    continue
  fi

  r1="${r1_candidates[0]}"
  r2="${r2_candidates[0]}"

  clean_r1="${OUTPUT_DIR}/${sample}_R1.clean.fq.gz"
  clean_r2="${OUTPUT_DIR}/${sample}_R2.clean.fq.gz"
  fastp_html="${OUTPUT_DIR}/${sample}_fastp.html"
  fastp_json="${OUTPUT_DIR}/${sample}_fastp.json"
  fastp_log="${OUTPUT_DIR}/${sample}_fastp.log"
  hisat2_log="${OUTPUT_DIR}/${sample}.hisat2.log"
  sorted_bam="${OUTPUT_DIR}/${sample}.sorted.bam"
  flagstat_file="${OUTPUT_DIR}/${sample}.flagstat.txt"

  echo "  fastp input R1: $r1"
  echo "  fastp input R2: $r2"

  fastp \
    --in1 "$r1" \
    --in2 "$r2" \
    --out1 "$clean_r1" \
    --out2 "$clean_r2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred "$QUALIFIED_QUALITY_PHRED" \
    --length_required "$LENGTH_REQUIRED" \
    --trim_poly_g \
    --html "$fastp_html" \
    --json "$fastp_json" \
    --thread "$FASTP_THREADS" \
    2> "$fastp_log"

  echo "  Aligning reads and creating a coordinate-sorted BAM file"

  hisat2 \
    -p "$HISAT2_THREADS" \
    -x "$HISAT2_INDEX" \
    -1 "$clean_r1" \
    -2 "$clean_r2" \
    --dta \
    --known-splicesite-infile "$SPLICE_SITES" \
    2> "$hisat2_log" \
  | samtools sort \
      -@ "$SAMTOOLS_THREADS" \
      -m "$SAMTOOLS_SORT_MEMORY" \
      -o "$sorted_bam" \
      -

  samtools index -@ "$SAMTOOLS_THREADS" "$sorted_bam"
  samtools flagstat -@ "$SAMTOOLS_THREADS" "$sorted_bam" > "$flagstat_file"

  ((processed_samples += 1))
  echo "Finished sample: $sample"
  echo
done

{
  echo "Analysis finished: $(date --iso-8601=seconds)"
  echo "Processed samples: $processed_samples"
  echo "Skipped samples: $skipped_samples"
} | tee "${OUTPUT_DIR}/pipeline_summary.txt"

echo "All available samples processed."
