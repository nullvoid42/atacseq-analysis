#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Pipeline Configuration
# Source this file: source config/config.sh
# ============================================================================

# Project settings
export PROJECT_NAME="GM12878_ATACseq"
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Reference genome
export GENOME="hg38"
export BOWTIE2_INDEX="/opt/genomes/hg38"
export CHROM_SIZES="/opt/genomes/hg38.chrom.sizes"

# Data directories
export DATA_DIR="$PROJECT_ROOT/data"
export RAW_DATA="$PROJECT_ROOT/results/raw_data"
export FASTQC_DIR="$PROJECT_ROOT/results/fastqc"
export TRIM_DIR="$PROJECT_ROOT/results/trimming"
export ALIGN_DIR="$PROJECT_ROOT/results/alignment"
export FILTER_DIR="$PROJECT_ROOT/results/alignment/filtered"
export PEAK_DIR="$PROJECT_ROOT/results/peaks"
export FOOTPRINT_DIR="$PROJECT_ROOT/results/footprints"
export PLOTS_DIR="$PROJECT_ROOT/results/plots"
export LOG_DIR="$PROJECT_ROOT/logs"

# Analysis parameters
export THREADS=8
export MAPQ_MIN=30
export GENOME_SIZE="hs"  # 3.0e9 for human
export EXT_SIZE=200
export QVALUE=0.01

# Adapter sequences (Nextera XT)
export ADAPTER_FWD="CTGTCTCTTATACACATCT"
export ADAPTER_REV="CTGTCTCTTATACACATCT"

# Sample configuration
export SAMPLE_ID="GM12878_ATAC_Rep1"
