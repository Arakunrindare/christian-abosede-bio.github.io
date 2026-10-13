#!/bin/bash
# =============================================================================
# rmats_splicing_analysis.sh
# Differential Alternative Splicing Analysis with rMATS
#
# Detects differential splicing events between two conditions from RNA-seq
# BAM files. Identifies five event types:
#   SE   — skipped exon
#   A5SS — alternative 5' splice site
#   A3SS — alternative 3' splice site
#   MXE  — mutually exclusive exons
#   RI   — retained intron
#
# Input:
#   Two plain-text files, each listing BAM paths comma-separated on one line:
#     GROUP1_BAMS — treatment / condition of interest
#     GROUP2_BAMS — control / reference condition
#
#   Generate these files from your BAM directory:
#     ls /path/to/mapped/*.bam | paste -sd, - > group1.txt
#     ls /path/to/mapped/*.bam | paste -sd, - > group2.txt
#
#   Check read length from a FASTQ file before running:
#     zcat sample_R1.fastq.gz | head -n 2 | tail -n 1 | wc -c
#   Subtract 1 (wc -c counts the newline) and set READ_LENGTH accordingly.
#
# Output:
#   OUTPUT_DIR/   — rMATS results (*.MATS.JC.txt, *.MATS.JCEC.txt per event)
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Run: bash rmats_splicing_analysis.sh
#
# Dependencies: rMATS (rmats.py), Python 3
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
GROUP1_BAMS=group1.txt              # File listing Group 1 BAM paths (comma-separated)
GROUP2_BAMS=group2.txt              # File listing Group 2 BAM paths (comma-separated)
GTF=/path/to/annotation.gtf        # Gene annotation GTF matching your genome
READ_LENGTH=90                      # Read length in bp (check with command above)
LIB_TYPE=fr-unstranded              # Library strandedness:
                                    #   fr-unstranded  — unstranded
                                    #   fr-firststrand — RF stranded (dUTP, SMARTer)
                                    #   fr-secondstrand — FR stranded
OUTPUT_DIR=rmats_output
THREADS=16
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/tmp_rmats"

rmats.py \
    --b1 "${GROUP1_BAMS}" \
    --b2 "${GROUP2_BAMS}" \
    --gtf "${GTF}" \
    --od "${OUTPUT_DIR}" \
    --tmp "${OUTPUT_DIR}/tmp_rmats" \
    --readLength "${READ_LENGTH}" \
    --variable-read-length \
    --libType "${LIB_TYPE}" \
    --nthread "${THREADS}"

echo "rMATS complete. Results in: ${OUTPUT_DIR}/"
echo "Next step: visualise significant events with rmats2sashimiplot.sh"
