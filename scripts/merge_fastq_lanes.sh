#!/bin/bash
# =============================================================================
# merge_fastq_lanes.sh
# Merge Paired-End FASTQ Files Across Sequencing Lanes or Batches
#
# Use case:
#   Samples sequenced across two lanes (L001, L002) or two sequencing runs
#   need their FASTQ files concatenated before alignment.
#   This script matches samples by name and merges R1/R2 files from both
#   sources into a single output pair per sample.
#
# Input directory structure:
#   BATCH1/   — contains L002 files: *_L002_R1_001.fastq.gz, *_L002_R2_001.fastq.gz
#   BATCH2/   — contains L001 files: *_L001_R1_001.fastq.gz, *_L001_R2_001.fastq.gz
#
# Output:
#   MERGED/   — one merged *_merged_R1.fastq.gz and *_merged_R2.fastq.gz per sample
#
# Usage:
#   1. Set BATCH1, BATCH2, and MERGED paths below.
#   2. Run: bash merge_fastq_lanes.sh
#
# Note: The order of concatenation (BATCH2 then BATCH1) places the earlier
#       lane first. Adjust if your batches have a different chronological order.
#
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
BATCH1="fastq-b1"   # Lane/batch containing L002 files
BATCH2="fastq-b2"   # Lane/batch containing L001 files
MERGED="merged"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${MERGED}"

for f1 in "${BATCH1}"/*_L002_R1_001.fastq.gz; do

  # Extract sample ID (everything before _L002)
  sample_id=$(basename "${f1}" _L002_R1_001.fastq.gz)

  # Locate all four input files
  r1_b1=$(find "${BATCH1}" -name "${sample_id}"'*_L002_R1_001.fastq.gz')
  r2_b1=$(find "${BATCH1}" -name "${sample_id}"'*_L002_R2_001.fastq.gz')
  r1_b2=$(find "${BATCH2}" -name "${sample_id}"'*_L001_R1_001.fastq.gz')
  r2_b2=$(find "${BATCH2}" -name "${sample_id}"'*_L001_R2_001.fastq.gz')

  # Output files
  out_r1="${MERGED}/${sample_id}_merged_R1.fastq.gz"
  out_r2="${MERGED}/${sample_id}_merged_R2.fastq.gz"

  # Merge only if all four files exist
  if [[ -f "${r1_b1}" && -f "${r1_b2}" && -f "${r2_b1}" && -f "${r2_b2}" ]]; then
    echo "Merging: ${sample_id}"
    cat "${r1_b2}" "${r1_b1}" > "${out_r1}"
    cat "${r2_b2}" "${r2_b1}" > "${out_r2}"
  else
    echo "WARNING: Missing file(s) for ${sample_id} — skipping."
    [[ ! -f "${r1_b1}" ]] && echo "  Missing: ${BATCH1}/${sample_id}*_L002_R1*"
    [[ ! -f "${r2_b1}" ]] && echo "  Missing: ${BATCH1}/${sample_id}*_L002_R2*"
    [[ ! -f "${r1_b2}" ]] && echo "  Missing: ${BATCH2}/${sample_id}*_L001_R1*"
    [[ ! -f "${r2_b2}" ]] && echo "  Missing: ${BATCH2}/${sample_id}*_L001_R2*"
  fi

done

echo "Merging complete. Output in: ${MERGED}/"
