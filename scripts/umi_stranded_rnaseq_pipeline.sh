#!/bin/bash
# =============================================================================
# umi_stranded_rnaseq_pipeline.sh
# UMI-Aware Stranded Total RNA-seq Pipeline
#
# Kit:          Takara SMARTer Stranded Total RNA-seq (Pico Input, Mammalian)
# Genome:       GRCm39 (mm39) — swap GENOME_DIR/GTF for any STAR-indexed genome
# Annotation:   GENCODE vM36 GTF
# UMI:          8 nt on Read 2 | Linker: 6 nt after UMI on Read 2
# Strandedness: RF (fr-firststrand) → featureCounts -s 2
#
# Pipeline steps:
#   1. QC of raw reads       (FastQC + MultiQC)
#   2. UMI extraction +
#      adapter/quality trim  (fastp)
#   3. QC after trimming     (FastQC + MultiQC)
#   4. Alignment             (STAR)
#   5. Alignment QC          (samtools flagstat + MultiQC)
#   6. UMI deduplication     (umi_tools dedup)
#   7. Read counting         (featureCounts)
#
# Usage:
#   1. Set the CONFIG variables below.
#   2. Place paired-end FASTQ files in RAW_DIR named *_R1_001.fastq.gz / *_R2_001.fastq.gz
#   3. Uncomment steps selectively or run all:
#        bash umi_stranded_rnaseq_pipeline.sh
#
# Dependencies: fastqc, multiqc, fastp, STAR, samtools, umi_tools, featureCounts
# Author: Christian Abosede (github.com/Arakunrindare)
# =============================================================================

# ── CONFIG ───────────────────────────────────────────────────────────────────
THREADS=16
GENOME_DIR=/path/to/STAR_genome_index          # STAR genome directory
GTF=/path/to/annotation.gtf                   # Gene annotation GTF
RAW_DIR=rawfiles                               # Directory with raw FASTQ files
TRIM_DIR=trimmed_fastp
MAPPED_DIR=mapped_STAR
DEDUP_DIR=dedup
COUNTS_DIR=featureCounts
# ─────────────────────────────────────────────────────────────────────────────


# =============================================================================
# STEP 1 — QUALITY CHECK (raw reads)
# =============================================================================

mkdir -p fastqc_raw

fastqc \
  -t ${THREADS} \
  -o fastqc_raw \
  ${RAW_DIR}/*.fastq.gz

multiqc \
  -o fastqc_raw \
  fastqc_raw/*.zip


# =============================================================================
# STEP 2 — UMI EXTRACTION + ADAPTER & QUALITY TRIMMING (fastp)
#
# Read 2 structure: [8 nt UMI][6 nt linker][insert]
#   --umi --umi_loc read2 --umi_len 8  → extracts UMI into read name
#   --trim_front2 14                   → removes 8 nt UMI + 6 nt linker from
#                                        R2 5' end (fastp retains UMI in the
#                                        read name after tagging, so trim all 14)
#   --adapter_sequence / _r2           → explicit TruSeq adapters (safer than
#                                        auto-detect alone for this kit)
#   --detect_adapter_for_pe            → overlap-based detection as a backup
#   --cut_right                        → sliding-window quality trimming
#   -l 30                              → discard reads shorter than 30 nt
# =============================================================================

mkdir -p ${TRIM_DIR}

for R1 in ${RAW_DIR}/*_R1_001.fastq.gz; do

  base=$(basename ${R1} _R1_001.fastq.gz)
  R2=${RAW_DIR}/${base}_R2_001.fastq.gz

  echo "=== Trimming: ${base} ==="

  fastp \
    -i ${R1} \
    -I ${R2} \
    -o ${TRIM_DIR}/${base}_R1.fastq.gz \
    -O ${TRIM_DIR}/${base}_R2.fastq.gz \
    --umi \
    --umi_loc read2 \
    --umi_len 8 \
    --trim_front2 14 \
    --adapter_sequence  AGATCGGAAGAGCACACGTCTGAAC \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGG \
    --detect_adapter_for_pe \
    --cut_right \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    -l 30 \
    -w ${THREADS} \
    -j ${TRIM_DIR}/${base}.fastp.json \
    -h ${TRIM_DIR}/${base}.fastp.html

done

# QC after trimming
mkdir -p ${TRIM_DIR}/fastqc

fastqc \
  -t ${THREADS} \
  -o ${TRIM_DIR}/fastqc \
  ${TRIM_DIR}/*.fastq.gz

multiqc \
  -o ${TRIM_DIR} \
  ${TRIM_DIR}/fastqc/*.zip


# =============================================================================
# STEP 3 — ALIGNMENT (STAR)
#
#   --outSAMattributes All              → includes NH, HI, nM needed by
#                                         umi_tools dedup
#   --outFilterMultimapNmax 20          → allow up to 20 multi-mappers
#   --outFilterMismatchNoverReadLmax 0.04 → 4% mismatch (Illumina standard)
#   --outFilterMatchNmin 25             → minimum 25 nt matching bases
#   --quantMode GeneCounts              → per-gene counts as QC sanity check
#                                         (column 4 = RF stranded)
#                                         Not used for final counts — use
#                                         featureCounts on deduplicated BAMs
# =============================================================================

mkdir -p ${MAPPED_DIR}

for R1 in ${TRIM_DIR}/*_R1.fastq.gz; do

  base=$(basename ${R1} _R1.fastq.gz)
  R2=${TRIM_DIR}/${base}_R2.fastq.gz

  echo "=== Mapping: ${base} ==="

  STAR \
    --genomeDir ${GENOME_DIR} \
    --runThreadN ${THREADS} \
    --readFilesIn ${R1} ${R2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${MAPPED_DIR}/${base}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes All \
    --quantMode GeneCounts \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 25

  samtools index \
    -@ ${THREADS} \
    ${MAPPED_DIR}/${base}_Aligned.sortedByCoord.out.bam

done

# Alignment QC
mkdir -p ${MAPPED_DIR}/flagstat

for BAM in ${MAPPED_DIR}/*_Aligned.sortedByCoord.out.bam; do
  base=$(basename ${BAM} _Aligned.sortedByCoord.out.bam)
  samtools flagstat \
    -@ ${THREADS} \
    ${BAM} > ${MAPPED_DIR}/flagstat/${base}.flagstat.txt
done

multiqc \
  -o ${MAPPED_DIR} \
  ${MAPPED_DIR}/


# =============================================================================
# STEP 4 — UMI DEDUPLICATION (umi_tools dedup)
#
# Must run BEFORE counting to remove PCR duplicates identified by UMI.
# fastp tagged UMIs into the read name during Step 2.
#   --paired              → paired-end mode
#   --unpaired-reads discard → discard reads whose mate was lost
# =============================================================================

mkdir -p ${DEDUP_DIR}

for BAM in ${MAPPED_DIR}/*_Aligned.sortedByCoord.out.bam; do

  base=$(basename ${BAM} _Aligned.sortedByCoord.out.bam)

  echo "=== Deduplicating: ${base} ==="

  umi_tools dedup \
    --paired \
    --unpaired-reads discard \
    -I ${BAM} \
    -S ${DEDUP_DIR}/${base}_dedup.bam \
    --log ${DEDUP_DIR}/${base}_dedup.log \
    --output-stats ${DEDUP_DIR}/${base}_dedup_stats

  samtools index \
    -@ ${THREADS} \
    ${DEDUP_DIR}/${base}_dedup.bam

done


# =============================================================================
# STEP 5 — READ COUNTING (featureCounts)
#
#   -s 2         → reversely stranded (RF / fr-firststrand) — REQUIRED for kit
#   -p           → paired-end mode
#   -t exon      → count at exon level
#   -g gene_name → summarise by gene name (change to gene_id if preferred)
#   -O           → allow reads overlapping multiple features
#   -M           → count multi-mapping reads
#   --fraction   → distribute multi-overlap/multi-mapper counts fractionally
#                  (use with -O and -M to avoid count inflation)
# Input: deduplicated BAMs from Step 4
# =============================================================================

mkdir -p ${COUNTS_DIR}

featureCounts \
  -p \
  -s 2 \
  -t exon \
  -g gene_name \
  -O \
  -M \
  --fraction \
  -T ${THREADS} \
  -a ${GTF} \
  -o ${COUNTS_DIR}/counts.txt \
  ${DEDUP_DIR}/*_dedup.bam

# Summary QC
multiqc \
  -o ${COUNTS_DIR} \
  ${COUNTS_DIR}/


# =============================================================================
# DONE
# Next steps (in R / DESeq2):
#   - Load counts.txt into DESeq2
#   - Filter baseMean < 10 before GO analysis
#   - Use adjusted p-value < 0.05 as significance threshold
#   - STAR GeneCounts column 4 can be used as sanity check
#     (RF stranded = column 4 of ReadsPerGene.out.tab)
# =============================================================================

echo "Pipeline complete."
