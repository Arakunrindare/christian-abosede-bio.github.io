#!/bin/bash
# =============================================================================
# bedgraph_to_bigwig.sh
# Convert bedGraph Files to BigWig Format for Genome Browser Visualisation
#
# Use case:
#   After producing genome coverage files in bedGraph format (e.g. from
#   bismark, bedtools genomecov, or deepTools), convert them to BigWig for
#   efficient loading in genome browsers (IGV, UCSC, etc.).
#
# Pipeline:
#   1. Decompress .bedGraph.gz files
#   2. Sort bedGraph by chromosome and start position
#   3. Convert sorted bedGraph to BigWig using chromosome sizes
#
# Input:
#   Directory of .bedGraph.gz files
#   Chromosome sizes file (chrom.sizes) matching your reference genome
#     → Generate with: fetchChromSizes <assembly> > genome.chrom.sizes
#       or: samtools faidx genome.fa && cut -f1,2 genome.fa.fai > genome.chrom.sizes
#
# Output:
#   bigwigs/   — one .bw file per input bedGraph
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Run: bash bedgraph_to_bigwig.sh
#
# Dependencies: bedGraphToBigWig (UCSC kent tools), sort (coreutils)
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
THREADS=8
INPUT_DIR=/path/to/bedgraph/files         # Directory containing .bedGraph.gz files
CHROM_SIZES=/path/to/genome.chrom.sizes   # Chromosome sizes file for your genome
SORTED_DIR=sorted_bedgraph
OUTPUT_DIR=bigwigs
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "${SORTED_DIR}" "${OUTPUT_DIR}"

for gz in "${INPUT_DIR}"/*.bedGraph.gz; do

  # Derive base name without extension
  base=$(basename "${gz}" .bedGraph.gz)

  echo "=== Processing: ${base} ==="

  # Step 1: Decompress
  gunzip -k "${gz}"   # -k keeps the original .gz file

  # Step 2: Sort by chromosome (lexicographic) then start position (numeric)
  sort -k1,1 -k2,2n \
    "${INPUT_DIR}/${base}.bedGraph" \
    > "${SORTED_DIR}/${base}_sorted.bedGraph"

  # Step 3: Convert to BigWig
  bedGraphToBigWig \
    "${SORTED_DIR}/${base}_sorted.bedGraph" \
    "${CHROM_SIZES}" \
    "${OUTPUT_DIR}/${base}.bw"

  # Clean up uncompressed bedGraph to save space
  rm -f "${INPUT_DIR}/${base}.bedGraph"

  echo "    Done: ${OUTPUT_DIR}/${base}.bw"

done

echo "All conversions complete. BigWig files are in: ${OUTPUT_DIR}/"
