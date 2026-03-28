#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 4: Filtering
# ============================================================================
# - Remove duplicates with Picard MarkDuplicates
# - Filter by MAPQ >= 30
# - Remove chrM (mitochondrial reads)
# - Proper pairs only
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
ALIGN_DIR="$PROJECT_ROOT/results/alignment"
FILTER_DIR="$PROJECT_ROOT/results/alignment/filtered"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$FILTER_DIR" "$LOG_DIR"

# Default parameters
MAPQ_MIN="${MAPQ_MIN:-30}"
REMOVE_CHRM="${REMOVE_CHRM:-true}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

remove_duplicates() {
    local input_bam="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Removing duplicates from: $sample"
    
    # Check for Picard
    if ! command -v picard &> /dev/null && ! command -v java &> /dev/null; then
        log "WARNING: Picard not found. Install with: conda install -c bioconda picard"
        log "Skipping duplicate removal..."
        echo "$input_bam"
        return 0
    fi
    
    local output_bam="$output_dir/${sample}.dedup.bam"
    local metrics="$output_dir/${sample}.dup_metrics.txt"
    
    # Try picard command first, then fall back to java -jar
    if command -v picard &> /dev/null; then
        picard MarkDuplicates \
            I="$input_bam" \
            O="$output_bam" \
            M="$metrics" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tee -a "$LOG_DIR/picard.log"
    else
        # Fallback to java directly
        java -jar "$PICARD_JAR" MarkDuplicates \
            I="$input_bam" \
            O="$output_bam" \
            M="$metrics" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tee -a "$LOG_DIR/picard.log"
    fi
    
    # Index
    samtools index "$output_bam"
    
    log "Duplicate removal complete: $output_bam"
    echo "$output_bam"
}

filter_bam() {
    local input_bam="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Filtering BAM: MAPQ >= $MAPQ_MIN, proper pairs, remove chrM"
    
    if ! command -v samtools &> /dev/null; then
        log "ERROR: SAMtools not found"
        return 1
    fi
    
    local output_bam="$output_dir/${sample}.filtered.bam"
    
    # Build samtools filter command
    # -F 4: keep mapped reads
    # -q $MAPQ_MIN: minimum mapping quality
    # -f 2: properly paired
    # -F 256: not secondary alignment
    
    local filter_flags="-F 4 -F 256 -f 2 -q $MAPQ_MIN"
    
    if [[ "$REMOVE_CHRM" == "true" ]]; then
        # Remove chrM and chrMT
        samtools view -b $filter_flags "$input_bam" \
            --exclude chrM \
            --exclude chrMT \
            -o "$output_bam" 2>&1 | tee -a "$LOG_DIR/filter.log"
    else
        samtools view -b $filter_flags "$input_bam" \
            -o "$output_bam" 2>&1 | tee -a "$LOG_DIR/filter.log"
    fi
    
    # Sort and index
    samtools sort "$output_bam" -o "${output_bam%.bam}.sorted.bam"
    samtools index "${output_bam%.bam}.sorted.bam"
    
    # Clean up unsorted
    rm -f "$output_bam" 2>/dev/null || true
    
    local final_bam="${output_bam%.bam}.sorted.bam"
    log "Filtering complete: $final_bam"
    echo "$final_bam"
}

fragment_size_distribution() {
    local bam_file="$1"
    local sample="$2"
    local output_file="$3"
    
    log "Calculating fragment size distribution..."
    
    if ! command -v samtools &> /dev/null; then
        log "WARNING: SAMtools not found"
        return 1
    fi
    
    # Get insert sizes from properly paired reads
    samtools view -f 2 "$bam_file" | \
        awk '{if($9>0) print $9}' | \
        sort | \
        uniq -c | \
        sort -b -k2,2n > "$output_file"
    
    log "Fragment sizes saved: $output_file"
}

generate_filter_report() {
    local sample="$1"
    local input_bam="$2"
    local dedup_bam="$3"
    local final_bam="$4"
    local output_file="$5"
    
    {
        echo "ATAC-seq Filtering Report: $sample"
        echo "=================================="
        echo "Date: $(date)"
        echo ""
        
        echo "Step 1: Duplicate Removal"
        echo "------------------------"
        echo "Input:  $input_bam"
        echo "Reads in:  $(samtools view -c "$input_bam" 2>/dev/null)"
        echo "Reads out: $(samtools view -c "$dedup_bam" 2>/dev/null)"
        echo ""
        
        echo "Step 2: Quality & Pair Filtering"
        echo "---------------------------------"
        echo "MAPQ threshold: >= $MAPQ_MIN"
        echo "Proper pairs: required"
        echo "chrM removal: $REMOVE_CHRM"
        echo "Input:  $dedup_bam"
        echo "Reads in:  $(samtools view -c "$dedup_bam" 2>/dev/null)"
        echo "Reads out: $(samtools view -c "$final_bam" 2>/dev/null)"
        echo ""
        
        echo "Final Statistics"
        echo "----------------"
        echo "Total filtered reads: $(samtools view -c "$final_bam" 2>/dev/null)"
        echo "Mapped (MAPQ>=30): $(samtools view -c -q $MAPQ_MIN "$final_bam" 2>/dev/null)"
        echo "Proper pairs: $(samtools view -c -f 2 "$final_bam" 2>/dev/null)"
        
    } > "$output_file"
    
    log "Filter report saved: $output_file"
}

# Main execution
main() {
    local sample_id="${1:-GM12878_ATAC_Rep1}"
    
    log "=== Step 4: Filtering ==="
    
    # Find input BAM from alignment step
    local input_bam="$ALIGN_DIR/${sample_id}.sorted.bam"
    
    if [[ ! -f "$input_bam" ]]; then
        log "ERROR: Input BAM not found: $input_bam"
        log "Please run alignment step first"
        exit 1
    fi
    
    log "Input BAM: $input_bam"
    
    # Step 4a: Remove duplicates
    local dedup_bam=$(remove_duplicates "$input_bam" "$sample_id" "$FILTER_DIR")
    
    # Step 4b: Filter by MAPQ and remove chrM
    local final_bam=$(filter_bam "$dedup_bam" "$sample_id" "$FILTER_DIR")
    
    # Step 4c: Fragment size distribution
    fragment_size_distribution "$final_bam" "$sample_id" \
        "$FILTER_DIR/${sample_id}.fragment_sizes.txt"
    
    # Step 4d: Generate report
    generate_filter_report "$sample_id" "$input_bam" "$dedup_bam" "$final_bam" \
        "$FILTER_DIR/${sample_id}_filter_report.txt"
    
    log "=== Step 4 Complete ==="
    log "Filtered BAM: $final_bam"
    log "Reports: $FILTER_DIR"
}

main "$@"
