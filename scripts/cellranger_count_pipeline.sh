#!/bin/bash
# =============================================================================
# cellranger_count_pipeline.sh
# Single-Cell RNA-seq Gene Expression Quantification with Cell Ranger
#
# Use case:
#   Process 10x Genomics Chromium single-cell RNA-seq (3' or 5' GEX) data
#   through the Cell Ranger count pipeline to produce per-cell gene-barcode
#   matrices for downstream analysis (Seurat, Scanpy, etc.).
#
# What Cell Ranger count does:
#   1. Aligns reads to the reference transcriptome (STAR internally)
#   2. Performs barcode correction and UMI deduplication
#   3. Generates a filtered feature-barcode matrix (cells × genes)
#   4. Produces a summary HTML report (web_summary.html)
#
# Output (inside ${OUTPUT_ID}/outs/):
#   filtered_feature_bc_matrix/   → use this for Seurat/Scanpy
#   raw_feature_bc_matrix/
#   possorted_genome_bam.bam
#   web_summary.html
#   metrics_summary.csv
#
# Usage:
#   1. Set the CONFIG variables below for each sample.
#   2. Run: bash cellranger_count_pipeline.sh
#   3. For multiple samples, copy the cellranger count block into the loop.
#
# Dependencies: Cell Ranger (≥7.0 recommended)
#               Reference transcriptome (pre-built or custom)
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
TRANSCRIPTOME=/path/to/cellranger/reference  # Pre-built Cell Ranger reference
                                              # Download from 10x Genomics or
                                              # build with cellranger mkref
FASTQ_DIR=/path/to/fastqs                    # Directory containing FASTQ files
SAMPLE_NAME=sample_name                      # Must match the sample prefix in
                                              # the FASTQ filenames (before _S*_L*)
OUTPUT_ID=cellranger_${SAMPLE_NAME}          # Unique run ID (output folder name)
# ─────────────────────────────────────────────────────────────────────────────

# Single-sample run
cellranger count \
  --id="${OUTPUT_ID}" \
  --transcriptome="${TRANSCRIPTOME}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE_NAME}" \
  --chemistry=auto    # auto-detect: threeprime, fiveprime, SC3Pv3, etc.
                      # Specify explicitly if auto-detection fails


# =============================================================================
# Multi-sample loop — uncomment and edit SAMPLES array for batch processing
# =============================================================================

# SAMPLES=("sample_A" "sample_B" "sample_C")
#
# for SAMPLE_NAME in "${SAMPLES[@]}"; do
#   OUTPUT_ID="cellranger_${SAMPLE_NAME}"
#   echo "=== Running Cell Ranger count: ${SAMPLE_NAME} ==="
#   cellranger count \
#     --id="${OUTPUT_ID}" \
#     --transcriptome="${TRANSCRIPTOME}" \
#     --fastqs="${FASTQ_DIR}" \
#     --sample="${SAMPLE_NAME}" \
#     --chemistry=auto
# done

echo "Cell Ranger count complete. Check ${OUTPUT_ID}/outs/web_summary.html"
