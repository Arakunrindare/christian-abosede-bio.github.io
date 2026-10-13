#!/bin/bash
# =============================================================================
# paired_end_rnaseq_pipeline.sh
# Standard Paired-End RNA-seq Pipeline
#
# Genome:  Any STAR-indexed genome (default: hg38)
# Trimmer: Trimmomatic PE
# Aligner: STAR
# Counter: featureCounts
#
# Pipeline steps:
#   1. QC of raw reads    (FastQC + MultiQC)
#   2. Adapter trimming   (Trimmomatic PE)
#   3. Genome indexing    (STAR — run once per genome)
#   4. Alignment          (STAR)
#   5. Read counting      (featureCounts)
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Create a plain-text file listing sample base names (one per line),
#      e.g. "sample1" for files sample1.R1.fastq / sample1.R2.fastq
#      and set SAMPLE_LIST to that file path.
#   3. Run: bash paired_end_rnaseq_pipeline.sh
#
# Dependencies: fastqc, multiqc, trimmomatic, STAR, samtools, featureCounts
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
THREADS=16
RAW_DIR=/path/to/raw/fastq          # Directory containing raw FASTQ files
GENOME_FASTA=/path/to/genome.fa     # Reference genome FASTA
GTF=/path/to/annotation.gtf        # Gene annotation GTF
ADAPTER_FA=TruSeq3-PE-2.fa          # Trimmomatic adapter file (in its share dir)
SAMPLE_LIST=samples.txt             # Plain-text file: one sample base name per line

TRIM_DIR=trimmed_fastq
GENOME_INDEX=genomeIndex
MAPPED_DIR=mapped
COUNTS_DIR=featureCounts
# ─────────────────────────────────────────────────────────────────────────────


# =============================================================================
# STEP 1 — QC OF RAW READS
# =============================================================================

mkdir -p fastqc_raw

find ${RAW_DIR} -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc \
  -t ${THREADS} \
  -o fastqc_raw

multiqc fastqc_raw -o fastqc_raw


# =============================================================================
# STEP 2 — ADAPTER & QUALITY TRIMMING (Trimmomatic PE)
#
#   ILLUMINACLIP:2:30:10  → remove TruSeq adapters (seed mismatches:2,
#                            palindrome clip threshold:30, simple clip:10)
#   HEADCROP:12           → remove first 12 bases (library-prep bias)
#   LEADING:3             → remove leading low-quality bases (< Q3)
#   TRAILING:3            → remove trailing low-quality bases (< Q3)
#   SLIDINGWINDOW:4:15    → sliding-window trimming (4 nt window, avg Q15)
#   MINLEN:30             → discard reads shorter than 30 nt after trimming
# =============================================================================

mkdir -p ${TRIM_DIR}

while IFS= read -r base; do
  echo "=== Trimming: ${base} ==="
  trimmomatic PE \
    -threads ${THREADS} \
    -phred33 \
    ${RAW_DIR}/${base}.R1.fastq \
    ${RAW_DIR}/${base}.R2.fastq \
    ${TRIM_DIR}/${base}.R1_trimmed_clean.fastq \
    ${TRIM_DIR}/${base}.R1_trimmed_unpaired.fastq \
    ${TRIM_DIR}/${base}.R2_trimmed_clean.fastq \
    ${TRIM_DIR}/${base}.R2_trimmed_unpaired.fastq \
    ILLUMINACLIP:${ADAPTER_FA}:2:30:10 \
    HEADCROP:12 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:30
done < ${SAMPLE_LIST}


# =============================================================================
# STEP 3 — BUILD GENOME INDEX (STAR)
# Run once per genome — skip if index already exists.
# Estimated time: ~30 min for hg38 with 16 threads.
# =============================================================================

mkdir -p ${GENOME_INDEX}

STAR \
  --runMode genomeGenerate \
  --genomeDir ${GENOME_INDEX} \
  --genomeFastaFiles ${GENOME_FASTA} \
  --sjdbGTFfile ${GTF} \
  --runThreadN ${THREADS}

wait


# =============================================================================
# STEP 4 — ALIGNMENT (STAR)
#
#   --outFilterMismatchNmax 3   → max 3 mismatches per read
#   --outFilterMultimapNmax 1   → keep uniquely mapping reads only
#   --outSAMattributes All      → include all optional BAM tags
# =============================================================================

mkdir -p ${MAPPED_DIR}

while IFS= read -r base; do
  echo "=== Mapping: ${base} ==="
  STAR \
    --genomeDir ${GENOME_INDEX} \
    --runThreadN ${THREADS} \
    --readFilesIn \
      ${TRIM_DIR}/${base}.R1_trimmed_clean.fastq \
      ${TRIM_DIR}/${base}.R2_trimmed_clean.fastq \
    --outFileNamePrefix ${MAPPED_DIR}/${base}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --outSAMattributes All \
    --outFilterMismatchNmax 3 \
    --outFilterMultimapNmax 1

  samtools index -@ ${THREADS} ${MAPPED_DIR}/${base}_Aligned.sortedByCoord.out.bam
done < ${SAMPLE_LIST}


# =============================================================================
# STEP 5 — READ COUNTING (featureCounts)
#
#   -O       → allow reads overlapping multiple features
#   -t gene  → count at gene level
# =============================================================================

mkdir -p ${COUNTS_DIR}

while IFS= read -r base; do
  echo "=== Counting: ${base} ==="
  featureCounts \
    -O \
    -t gene \
    -a ${GTF} \
    -o ${COUNTS_DIR}/counts_${base}.txt \
    ${MAPPED_DIR}/${base}*.bam
done < ${SAMPLE_LIST}

multiqc ${COUNTS_DIR} -o ${COUNTS_DIR}

echo "Pipeline complete."
