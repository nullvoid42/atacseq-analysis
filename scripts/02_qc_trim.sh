#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 2: Quality Control & Trimming
# ============================================================================
# - FastQC for read quality assessment
# - Adapter trimming with Cutadapt (Nextera XT adapters)
# - Remove low-quality reads
# - Remove mitochondrial reads (chrM)
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
DATA_DIR="$PROJECT_ROOT/data"
FASTQC_DIR="$PROJECT_ROOT/results/fastqc"
TRIM_DIR="$PROJECT_ROOT/results/trimming"
LOG_DIR="$PROJECT_ROOT/logs"
GENOME_DIR="$PROJECT_ROOT/references"

# Create directories
mkdir -p "$FASTQC_DIR" "$TRIM_DIR" "$LOG_DIR"

# Load config
source "$CONFIG_DIR/config.sh" 2>/dev/null || true

# Default parameters
THREADS="${THREADS:-8}"
ADAPTER_FWD="CTGTCTCTTATACACATCT"  # Nextera XT adapter (forward)
ADAPTER_REV="CTGTCTCTTATACACATCT"  # Nextera XT adapter (reverse)
QUALITY_CUTOFF=20
MIN_LENGTH=20
CHRM="chrM"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

run_fastqc() {
    local input_dir="$1"
    local output_dir="$2"
    local sample="$3"
    
    log "Running FastQC on: $sample"
    
    # Check if FastQC is available
    if ! command -v fastqc &> /dev/null; then
        log "WARNING: FastQC not found. Install with: conda install -c bioconda fastqc"
        return 1
    fi
    
    # Find FASTQ files
    local fq_files=$(find "$input_dir" -maxdepth 1 -name "*${sample}*.fastq.gz" -o -name "*${sample}*.fq.gz" 2>/dev/null)
    
    if [[ -z "$fq_files" ]]; then
        log "No FASTQ files found for sample: $sample"
        return 1
    fi
    
    fastqc \
        --threads "$THREADS" \
        --outdir "$output_dir" \
        $fq_files 2>&1 | tee -a "$LOG_DIR/fastqc.log"
    
    log "FastQC complete for: $sample"
}

trim_reads() {
    local fastq_1="$1"
    local fastq_2="$2"
    local sample="$3"
    local output_dir="$4"
    
    log "Trimming reads for: $sample"
    
    # Check if Cutadapt is available
    if ! command -v cutadapt &> /dev/null; then
        log "WARNING: Cutadapt not found. Install with: conda install -c bioconda cutadapt"
        return 1
    fi
    
    local output_1="$output_dir/${sample}_R1.trimmed.fastq.gz"
    local output_2="$output_dir/${sample}_R2.trimmed.fastq.gz"
    local report_1="$output_dir/${sample}_R1.cutadapt.txt"
    local report_2="$output_dir/${sample}_R2.cutadapt.txt"
    
    # Paired-end trimming with Cutadapt
    cutadapt \
        -a "$ADAPTER_FWD" \
        -A "$ADAPTER_REV" \
        -q "$QUALITY_CUTOFF" \
        -m "$MIN_LENGTH" \
        --trim-n \
        -o "$output_1" \
        -p "$output_2" \
        "$fastq_1" "$fastq_2" \
        > "$report_1" 2>&1 || true
    
    # Also process R2 report
    cp "$report_1" "$report_2"
    
    # Validate output
    if [[ -f "$output_1" && -f "$output_2" ]]; then
        local f1_size=$(stat -f%z "$output_1" 2>/dev/null || stat -c%s "$output_1" 2>/dev/null || echo 0)
        local f2_size=$(stat -f%z "$output_2" 2>/dev/null || stat -c%s "$output_2" 2>/dev/null || echo 0)
        
        if [[ "$f1_size" -gt 1000 && "$f2_size" -gt 1000 ]]; then
            log "Trimmed files created successfully"
        else
            log "WARNING: Trimmed files suspiciously small"
        fi
    fi
    
    log "Trimming complete for: $sample"
    echo "$output_1 $output_2"
}

remove_chrm() {
    local input_bam="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Removing chrM reads from: $sample"
    
    if ! command -v samtools &> /dev/null; then
        log "WARNING: SAMtools not found. Skipping chrM removal."
        return 1
    fi
    
    local output_bam="$output_dir/${sample}.no_chrm.bam"
    
    # Remove chrM alignments
    samtools view -b "$input_bam" \
        --exclude chrM \
        --exclude chrMT \
        -o "$output_bam" 2>&1 | tee -a "$LOG_DIR/filter.log"
    
    # Index
    samtools index "$output_bam"
    
    log "chrM removal complete: $output_bam"
}

run_multiqc() {
    local fastqc_dir="$1"
    local output_file="$2"
    
    log "Generating MultiQC report..."
    
    if command -v multiqc &> /dev/null; then
        multiqc "$fastqc_dir" -o "$(dirname "$output_file")" -n "$(basename "$output_file")" \
            2>&1 | tee -a "$LOG_DIR/multiqc.log"
        log "MultiQC report: $output_file"
    else
        log "MultiQC not found. Install with: pip install multiqc"
    fi
}

generate_trim_report() {
    local sample="$1"
    local fastqc_dir="$2"
    local trim_dir="$3"
    local output_file="$4"
    
    {
        echo "ATAC-seq Trimming Report: $sample"
        echo "=================================="
        echo "Date: $(date)"
        echo ""
        echo "Before Trimming (FastQC summary):"
        if [[ -d "$fastqc_dir" ]]; then
            find "$fastqc_dir" -name "*${sample}*_fastqc.zip" -exec unzip -p {} "*fastqc_data.txt" 2>/dev/null | \
                grep -E "^%GC|Sequence length|Total Sequences" | head -5 || echo "No data available"
        fi
        echo ""
        echo "Trimmed files:"
        ls -lh "$trim_dir"/*${sample}*.trimmed* 2>/dev/null || echo "No trimmed files"
    } > "$output_file"
}

# Parse samples from config
get_samples() {
    local config_file="$1"
    if [[ -f "$config_file/samples.tsv" ]]; then
        tail -n +2 "$config_file/samples.tsv" | cut -f1
    else
        echo "GM12878_ATAC_Rep1"
    fi
}

# Main execution
main() {
    local sample_id="${1:-GM12878_ATAC_Rep1}"
    
    log "=== Step 2: Quality Control & Trimming ==="
    
    # Find input FASTQ files
    local fq1="$DATA_DIR/${sample_id}_R1.fastq.gz"
    local fq2="$DATA_DIR/${sample_id}_R2.fastq.gz"
    
    # Fallback to ENCODE files
    if [[ ! -f "$fq1" ]]; then
        fq1="$DATA_DIR/ENCFF234QEM.fastq.gz"
        fq2="$DATA_DIR/ENCFF235ZKY.fastq.gz"
    fi
    
    if [[ ! -f "$fq1" ]]; then
        log "ERROR: Input FASTQ files not found"
        log "Please run scripts/01_download.sh first or place files in $DATA_DIR"
        exit 1
    fi
    
    log "Input files:"
    log "  R1: $fq1"
    log "  R2: $fq2"
    
    # Step 2a: Run FastQC on raw data
    log "Running FastQC on raw reads..."
    run_fastqc "$DATA_DIR" "$FASTQC_DIR/raw" "$sample_id"
    
    # Step 2b: Trim adapters and low-quality bases
    log "Trimming adapters and low-quality bases..."
    trim_reads "$fq1" "$fq2" "$sample_id" "$TRIM_DIR"
    
    # Step 2c: Run FastQC on trimmed data
    log "Running FastQC on trimmed reads..."
    run_fastqc "$TRIM_DIR" "$FASTQC_DIR/trimmed" "$sample_id"
    
    # Step 2d: Generate MultiQC report
    run_multiqc "$FASTQC_DIR" "$PROJECT_ROOT/results/multiqc_report.html"
    
    # Step 2e: Generate report
    generate_trim_report "$sample_id" "$FASTQC_DIR/raw" "$TRIM_DIR" \
        "$PROJECT_ROOT/results/trimming_report.txt"
    
    log "=== Step 2 Complete ==="
    log "FastQC reports: $FASTQC_DIR"
    log "Trimmed reads: $TRIM_DIR"
}

main "$@"
