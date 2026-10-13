#!/bin/bash
# =============================================================================
# rmats2sashimiplot_visualisation.sh
# Sashimi Plot Visualisation of rMATS Splicing Events
#
# Generates sashimi plots from rMATS output to visualise differential splicing
# events at specific genomic coordinates. Typically run after
# rmats_splicing_analysis.sh to visualise significant events of interest.
#
# Event types (set EVENT_TYPE to one of):
#   SE   — skipped exon          (most common; start here)
#   A5SS — alternative 5' splice site
#   A3SS — alternative 3' splice site
#   MXE  — mutually exclusive exons
#   RI   — retained intron
#
# Finding coordinates for EXON_S and INTRON_S:
#   Open the rMATS output file (e.g. SE.MATS.JC.txt) and find the row for
#   your event of interest. EXON_S is the exon scale (typically the exon
#   length in bp) and INTRON_S is the intron shrink scale (typically 5–20).
#   These control how the sashimi plot is scaled, not the genomic position.
#   Reasonable starting values: EXON_S=1  INTRON_S=5
#
# Group info file (GROUP_INFO):
#   Optional tab-separated file defining sample groupings for colour-coding.
#   Format (one sample per line):
#     SampleName <tab> GroupLabel
#   Leave GROUP_INFO empty ("") to skip grouping.
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Run: bash rmats2sashimiplot_visualisation.sh
#
# Dependencies: rmats2sashimiplot, Python 3, matplotlib, samtools
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
OUTPUT_DIR=sashimi_output              # Output directory for sashimi plots

LABEL_GROUP1="Condition_1"             # Display label for Group 1 (e.g. WT, Control)
LABEL_GROUP2="Condition_2"             # Display label for Group 2 (e.g. KO, Treatment)

EVENT_TYPE=SE                          # Splicing event type: SE A5SS A3SS MXE RI

RMATS_EVENT_FILE=/path/to/rmats_output/SE.MATS.JC.txt
                                       # Full path to the rMATS event file
                                       # Change SE to match your EVENT_TYPE above

GROUP1_BAMS=group1.txt                 # File listing Group 1 BAM paths (comma-separated)
GROUP2_BAMS=group2.txt                 # File listing Group 2 BAM paths (comma-separated)

EXON_S=1                               # Exon scaling factor (bp per unit; start with 1)
INTRON_S=5                             # Intron shrink scale (5–20; larger = more compressed)

GROUP_INFO=grouping.txt                # Sample grouping file (tab-separated; leave "" to skip)
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${OUTPUT_DIR}"

if [[ -n "${GROUP_INFO}" && -f "${GROUP_INFO}" ]]; then
    rmats2sashimiplot \
        -o "${OUTPUT_DIR}" \
        --l1 "${LABEL_GROUP1}" \
        --l2 "${LABEL_GROUP2}" \
        --event-type "${EVENT_TYPE}" \
        -e "${RMATS_EVENT_FILE}" \
        --b1 "${GROUP1_BAMS}" \
        --b2 "${GROUP2_BAMS}" \
        --exon_s "${EXON_S}" \
        --intron_s "${INTRON_S}" \
        --group-info "${GROUP_INFO}"
else
    rmats2sashimiplot \
        -o "${OUTPUT_DIR}" \
        --l1 "${LABEL_GROUP1}" \
        --l2 "${LABEL_GROUP2}" \
        --event-type "${EVENT_TYPE}" \
        -e "${RMATS_EVENT_FILE}" \
        --b1 "${GROUP1_BAMS}" \
        --b2 "${GROUP2_BAMS}" \
        --exon_s "${EXON_S}" \
        --intron_s "${INTRON_S}"
fi

echo "Sashimi plots complete. Output in: ${OUTPUT_DIR}/"
