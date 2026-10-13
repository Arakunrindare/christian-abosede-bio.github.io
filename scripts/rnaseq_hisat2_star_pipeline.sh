#!/bin/bash
# =============================================================================
# rnaseq_hisat2_star_pipeline.sh
# RNA-seq Pipeline with Choice of Aligner: HISAT2 or STAR
#
# Genome:  Any indexed genome (configure paths below)
# Trimmer: Trimmomatic PE or cutadapt (both shown — uncomment preferred)
# Aligner: HISAT2 (Section A) or STAR (Section B) — uncomment one
# Counter: featureCounts or HTSeq-count (both shown)
#
# This script demonstrates two common RNA-seq alignment strategies side by side.
# HISAT2 is splice-aware and fast; STAR handles longer reads and complex splicing.
# Both produce sorted BAMs compatible with featureCounts and HTSeq.
#
# Pipeline steps:
#   1.  Download public data from SRA   (fasterq-dump — optional)
#   2.  QC of raw reads                 (FastQC + MultiQC)
#   3A. Trimming with Trimmomatic PE    — or —
#   3B. Trimming with cutadapt
#   4.  QC after trimming               (FastQC + MultiQC)
#   5.  Build HISAT2 genome index       (run once)
#   6A. Alignment with HISAT2           — or —
#   6B. Alignment with STAR
#   7A. Read counting with featureCounts — or —
#   7B. Read counting with HTSeq-count
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Uncomment the trimmer and aligner you want.
#   3. Run: bash rnaseq_hisat2_star_pipeline.sh
#
# Dependencies: sra-tools, fastqc, multiqc, trimmomatic or cutadapt,
#               hisat2 or STAR, samtools, featureCounts, htseq-count
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
THREADS=64
GENOME_FASTA=/path/to/genome.fa         # Reference genome FASTA
GENOME_INDEX_HISAT2=/path/to/hisat2_index/genome_index
GENOME_INDEX_STAR=/path/to/STAR_genome_index
GTF=/path/to/annotation.gtf
SRR_LIST=SRR_Acc_List.txt               # SRA accession list (one per line)
RAW_DIR=rawdata_fastq
MERGED_DIR=merged
TRIM_TRIMMOMATIC=trimmed_trimmomatic
TRIM_CUTADAPT=trimmed_cutadapt
MAPPED_HISAT2=mapped_HISAT2
MAPPED_STAR=mapped_STAR
COUNTS_DIR=featureCounts
HTSEQ_DIR=HTseq_counts
# ─────────────────────────────────────────────────────────────────────────────


# =============================================================================
# STEP 1 — DOWNLOAD PUBLIC DATA FROM SRA (optional)
# Skip this step if you already have FASTQ files locally.
# =============================================================================

# mkdir -p ${RAW_DIR}
# while IFS= read -r line; do
#   echo "Downloading: ${line}"
#   fasterq-dump --outdir ${RAW_DIR} --split-3 ${line}
#   gzip ${RAW_DIR}/${line}*.fastq
# done < ${SRR_LIST}


# =============================================================================
# STEP 2 — QC OF RAW READS
# =============================================================================

# mkdir -p fastqc_raw
# fastqc -t ${THREADS} -o fastqc_raw ${MERGED_DIR}/*.fastq.gz
# multiqc fastqc_raw -o fastqc_raw


# =============================================================================
# STEP 3A — TRIMMING WITH TRIMMOMATIC PE
# Uncomment this block OR block 3B (cutadapt) — not both.
#
#   HEADCROP:12           → remove first 12 bases (library-prep bias)
#   LEADING/TRAILING:10   → remove low-quality ends (< Q10)
#   SLIDINGWINDOW:4:30    → sliding-window trimming (4 nt, avg Q30)
#   MINLEN:30             → discard reads shorter than 30 nt
# =============================================================================

# mkdir -p ${TRIM_TRIMMOMATIC}
# for R1 in ${MERGED_DIR}/*_R1.fastq.gz; do
#   base=$(basename ${R1} _R1.fastq.gz)
#   trimmomatic PE \
#     -threads ${THREADS} \
#     -phred33 \
#     ${MERGED_DIR}/${base}_R1.fastq.gz \
#     ${MERGED_DIR}/${base}_R2.fastq.gz \
#     ${TRIM_TRIMMOMATIC}/${base}_R1.clean.fastq.gz \
#     ${TRIM_TRIMMOMATIC}/${base}_R1.unpaired.fastq.gz \
#     ${TRIM_TRIMMOMATIC}/${base}_R2.clean.fastq.gz \
#     ${TRIM_TRIMMOMATIC}/${base}_R2.unpaired.fastq.gz \
#     ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
#     HEADCROP:12 \
#     LEADING:10 \
#     TRAILING:10 \
#     SLIDINGWINDOW:4:30 \
#     MINLEN:30
# done


# =============================================================================
# STEP 3B — TRIMMING WITH CUTADAPT (alternative to 3A)
#
#   --cut 10              → hard-clip first 10 bases from each read (5' bias)
#   --quality-cutoff 30,30 → trim 3' and 5' ends below Q30
#   --minimum-length 30   → discard reads shorter than 30 nt
# =============================================================================

# mkdir -p ${TRIM_CUTADAPT}
# for R1 in ${RAW_DIR}/*_R1.fastq.gz; do
#   base=$(basename ${R1} _R1.fastq.gz)
#   cutadapt \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     --cut 10 \
#     -j ${THREADS} \
#     --quality-cutoff 30,30 \
#     --minimum-length 30 \
#     -o ${TRIM_CUTADAPT}/${base}_R1.trimmed.fastq.gz \
#     -p ${TRIM_CUTADAPT}/${base}_R2.trimmed.fastq.gz \
#     ${RAW_DIR}/${base}_R1.fastq.gz \
#     ${RAW_DIR}/${base}_R2.fastq.gz
# done


# =============================================================================
# STEP 4 — QC AFTER TRIMMING
# =============================================================================

# mkdir -p ${TRIM_TRIMMOMATIC}/fastqc
# fastqc -t ${THREADS} -o ${TRIM_TRIMMOMATIC}/fastqc ${TRIM_TRIMMOMATIC}/*.clean.fastq.gz
# multiqc ${TRIM_TRIMMOMATIC}/fastqc -o ${TRIM_TRIMMOMATIC}


# =============================================================================
# STEP 5 — BUILD HISAT2 GENOME INDEX (run once per genome)
# =============================================================================

# hisat2-build ${GENOME_FASTA} ${GENOME_INDEX_HISAT2}


# =============================================================================
# STEP 6A — ALIGNMENT WITH HISAT2
# Uncomment this block OR block 6B (STAR) — not both.
# =============================================================================

# TRIM_DIR=${TRIM_TRIMMOMATIC}   # or ${TRIM_CUTADAPT}
# mkdir -p ${MAPPED_HISAT2}
# for R1 in ${TRIM_DIR}/*_R1.clean.fastq.gz; do
#   base=$(basename ${R1} _R1.clean.fastq.gz)
#   hisat2 \
#     -x ${GENOME_INDEX_HISAT2} \
#     -1 ${TRIM_DIR}/${base}_R1.clean.fastq.gz \
#     -2 ${TRIM_DIR}/${base}_R2.clean.fastq.gz \
#     -p ${THREADS} \
#     --dta \
#     | samtools sort -@ ${THREADS} -o ${MAPPED_HISAT2}/${base}.sorted.bam
#   samtools index -@ ${THREADS} ${MAPPED_HISAT2}/${base}.sorted.bam
# done


# =============================================================================
# STEP 6B — ALIGNMENT WITH STAR (alternative to 6A)
# =============================================================================

# TRIM_DIR=${TRIM_TRIMMOMATIC}   # or ${TRIM_CUTADAPT}
# mkdir -p ${MAPPED_STAR}
# for R1 in ${TRIM_DIR}/*_R1.clean.fastq.gz; do
#   base=$(basename ${R1} _R1.clean.fastq.gz)
#   STAR \
#     --genomeDir ${GENOME_INDEX_STAR} \
#     --runThreadN ${THREADS} \
#     --readFilesIn ${TRIM_DIR}/${base}_R1.clean.fastq.gz \
#                   ${TRIM_DIR}/${base}_R2.clean.fastq.gz \
#     --outFileNamePrefix ${MAPPED_STAR}/${base}_ \
#     --outSAMtype BAM SortedByCoordinate \
#     --outSAMunmapped None \
#     --outFilterMismatchNmax 3 \
#     --outFilterMultimapNmax 1 \
#     --outSAMattributes All \
#     --readFilesCommand zcat
#   samtools index -@ ${THREADS} ${MAPPED_STAR}/${base}_Aligned.sortedByCoord.out.bam
# done


# =============================================================================
# STEP 7A — READ COUNTING WITH featureCounts
#
#   -p      → paired-end mode
#   -t exon → count at exon level
#   -g gene_id → summarise by Ensembl gene ID
# =============================================================================

mkdir -p ${COUNTS_DIR}

# For HISAT2 BAMs:
# for BAM in ${MAPPED_HISAT2}/*.sorted.bam; do
#   base=$(basename ${BAM} .sorted.bam)
#   featureCounts -p -T ${THREADS} -t exon -g gene_id \
#     -a ${GTF} \
#     -o ${COUNTS_DIR}/counts_${base}.txt \
#     ${BAM}
# done

# For STAR BAMs:
featureCounts \
  -p \
  -T ${THREADS} \
  -t exon \
  -g gene_id \
  -a ${GTF} \
  -o ${COUNTS_DIR}/counts_STAR.txt \
  ${MAPPED_STAR}/*.bam


# =============================================================================
# STEP 7B — READ COUNTING WITH HTSeq-count (alternative to 7A)
#
#   -f bam    → input format
#   -r pos    → reads sorted by position
#   -s no     → unstranded (change to yes/reverse as appropriate)
#   -i gene_id → feature attribute to summarise by
# =============================================================================

# mkdir -p ${HTSEQ_DIR}
# for BAM in ${MAPPED_HISAT2}/*.sorted.bam; do
#   base=$(basename ${BAM} .sorted.bam)
#   htseq-count \
#     -f bam \
#     -r pos \
#     -s no \
#     -i gene_id \
#     ${BAM} \
#     ${GTF} > ${HTSEQ_DIR}/counts_${base}_HTseq.txt
# done

echo "Pipeline complete."
