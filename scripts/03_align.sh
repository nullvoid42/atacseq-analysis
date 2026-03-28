#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 3: Alignment
# ============================================================================
# - Align trimmed reads to hg38 with Bowtie2 (very-sensitive mode)
# - Convert to BAM, sort, and index using SAMtools
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
DATA_DIR="$PROJECT_ROOT/data"
TRIM_DIR="$PROJECT_ROOT/results/trimming"
ALIGN_DIR="$PROJECT_ROOT/results/alignment"
LOG_DIR="$PROJECT_ROOT/logs"
GENOME_INDEX="${BOWTIE2_INDEX:-/opt/genomes/hg38}"

mkdir -p "$ALIGN_DIR" "$LOG_DIR"

# Default parameters
THREADS="${THREADS:-8}"
MODE="very-sensitive"
MAX_FRAGMENT_LENGTH=2000
ORIENTATION="--fr"  # Forward-reverse for paired-end

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_bowtie2_index() {
    local index_base="$1"
    
    # Check for index files
    if [[ ! -f "${index_base}.1.bt2" ]]; then
        log "WARNING: Bowtie2 index not found at: $index_base"
        log "Please download hg38 index:"
        log "  bowtie2-build -f hg38.fa hg38"
        return 1
    fi
    return 0
}

align_with_bowtie2() {
    local fastq_1="$1"
    local fastq_2="$2"
    local sample="$3"
    local index_base="$4"
    local output_dir="$5"
    
    log "Aligning reads for: $sample"
    
    # Check for Bowtie2
    if ! command -v bowtie2 &> /dev/null; then
        log "ERROR: Bowtie2 not found. Install with: conda install -c bioconda bowtie2"
        return 1
    fi
    
    local sam_file="$output_dir/${sample}.aligned.sam"
    local log_file="$output_dir/${sample}.bowtie2.log"
    
    # Paired-end alignment with Bowtie2
    # very-sensitive mode: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    bowtie2 \
        --$MODE \
        -X "$MAX_FRAGMENT_LENGTH" \
        $ORIENTATION \
        -x "$index_base" \
        -1 "$fastq_1" \
        -2 "$fastq_2" \
        -S "$sam_file" \
        --threads "$THREADS" \
        --met-file "$log_file" \
        2>&1 | tee -a "$log_file"
    
    # Check alignment success
    if [[ ! -f "$sam_file" || $(wc -l < "$sam_file") -lt 10 ]]; then
        log "ERROR: Alignment failed or produced empty output"
        return 1
    fi
    
    log "Alignment complete: $sam_file"
    echo "$sam_file"
}

sam_to_bam() {
    local sam_file="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Converting SAM to BAM..."
    
    if ! command -v samtools &> /dev/null; then
        log "ERROR: SAMtools not found. Install with: conda install -c bioconda samtools"
        return 1
    fi
    
    local bam_file="$output_dir/${sample}.bam"
    
    # Convert SAM to BAM
    samtools view -bS "$sam_file" > "$bam_file" 2>/dev/null
    
    # Sort BAM
    samtools sort "$bam_file" -o "${bam_file%.bam}.sorted.bam"
    
    # Index
    samtools index "${bam_file%.bam}.sorted.bam"
    
    # Clean up unsorted BAM and SAM
    rm -f "$sam_file" "$bam_file" 2>/dev/null || true
    
    local sorted_bam="${bam_file%.bam}.sorted.bam"
    log "BAM created: $sorted_bam"
    echo "$sorted_bam"
}

alignment_stats() {
    local bam_file="$1"
    local sample="$2"
    local output_file="$3"
    
    log "Generating alignment statistics..."
    
    {
        echo "Alignment Statistics: $sample"
        echo "==========================="
        echo "Date: $(date)"
        echo ""
        echo "File: $bam_file"
        echo ""
        
        echo "Total alignments:"
        samtools view -c "$bam_file" 2>/dev/null | awk '{printf "  %d total reads\n", $1}'
        
        echo "Mapped reads:"
        samtools view -c -F 4 "$bam_file" 2>/dev/null | awk '{printf "  %d mapped\n", $1}'
        
        echo "Unmapped reads:"
        samtools view -c -f 4 "$bam_file" 2>/dev/null | awk '{printf "  %d unmapped\n", $1}'
        
        echo "Properly paired reads:"
        samtools view -c -f 2 -F 4 "$bam_file" 2>/dev/null | awk '{printf "  %d properly paired\n", $1}'
        
        echo ""
        echo "Alignment quality distribution:"
        samtools view "$bam_file" | awk '{print $5}' | sort | uniq -c | sort -rn | head -20
        
    } > "$output_file"
    
    log "Statistics saved: $output_file"
}

# Main execution
main() {
    local sample_id="${1:-GM12878_ATAC_Rep1}"
    
    log "=== Step 3: Alignment ==="
    
    # Find input FASTQ files (trimmed)
    local fq1="$TRIM_DIR/${sample_id}_R1.trimmed.fastq.gz"
    local fq2="$TRIM_DIR/${sample_id}_R2.trimmed.fastq.gz"
    
    # Fallback to raw data if trimmed not available
    if [[ ! -f "$fq1" ]]; then
        fq1="$DATA_DIR/${sample_id}_R1.fastq.gz"
        fq2="$DATA_DIR/${sample_id}_R2.fastq.gz"
    fi
    
    # Fallback to ENCODE files
    if [[ ! -f "$fq1" ]]; then
        fq1="$DATA_DIR/ENCFF234QEM.fastq.gz"
        fq2="$DATA_DIR/ENCFF235ZKY.fastq.gz"
    fi
    
    if [[ ! -f "$fq1" ]]; then
        log "ERROR: Input FASTQ files not found"
        log "Please run download and trimming steps first"
        exit 1
    fi
    
    log "Input files:"
    log "  R1: $fq1"
    log "  R2: $fq2"
    log "  Index: $GENOME_INDEX"
    
    # Check index
    check_bowtie2_index "$GENOME_INDEX" || true
    
    # Align
    local sam_file=$(align_with_bowtie2 "$fq1" "$fq2" "$sample_id" "$GENOME_INDEX" "$ALIGN_DIR")
    
    # Convert to sorted BAM
    local bam_file=$(sam_to_bam "$sam_file" "$sample_id" "$ALIGN_DIR")
    
    # Generate statistics
    alignment_stats "$bam_file" "$sample_id" "$ALIGN_DIR/${sample_id}_alignment_stats.txt"
    
    log "=== Step 3 Complete ==="
    log "Sorted BAM: $bam_file"
}

main "$@"
