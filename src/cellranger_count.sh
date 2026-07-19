#!/bin/bash
# Cell Ranger count pipeline for 10x scRNA-seq data (Rat mRatBN7.2)

# --- Download and prepare reference genome (run once) ---
# Reference: 10x Genomics pre-built reference for mRatBN7.2
# wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mRatBN7-2-2024-A.tar.gz"
# tar -xzvf refdata-gex-mRatBN7-2-2024-A.tar.gz
REFERENCE_DIR="/home/ubuntu/10x_genomics/references/refdata-gex-mRatBN7-2-2024-A"

# --- Define common parameters ---
LOCALCORES=16
LOCALMEM=64
EXPECT_CELLS=5000

# --- Process each sample ---
# Sample: BHS1
SAMPLE="BHS1"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/BHS1"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/BHS1"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=$LOCALMEM \
    --expect-cells=$EXPECT_CELLS

# Sample: BHS2
SAMPLE="BHS2"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/BHS2"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/BHS2"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=$LOCALMEM \
    --expect-cells=$EXPECT_CELLS

# Sample: BLS1
SAMPLE="BLS1"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/BLS1"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/BLS1"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS

# Sample: BLS2
SAMPLE="BLS2"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/BLS2"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/BLS2"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS

# Sample: KHS1
SAMPLE="KHS1"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/KHS1"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/KHS1"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS

# Sample: KHS2
SAMPLE="KHS2"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/KHS2"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/KHS2"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS

# Sample: KLS1
SAMPLE="KLS1"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/KLS1"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/KLS1"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS

# Sample: KLS2
SAMPLE="KLS2"
FASTQ_DIR="/home/ubuntu/DATA2/sc_RNA_seq/KLS2"
OUTPUT_DIR="/home/ubuntu/DATA2/sc_RNA_seq/processed_data/KLS2"
mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"
cellranger count \
    --id="$SAMPLE" \
    --create-bam=true \
    --transcriptome="$REFERENCE_DIR" \
    --fastqs="$FASTQ_DIR" \
    --sample="$SAMPLE" \
    --localcores=$LOCALCORES \
    --localmem=32 \
    --expect-cells=$EXPECT_CELLS