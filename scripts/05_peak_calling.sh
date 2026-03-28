#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 5: Peak Calling
# ============================================================================
# - Call peaks using MACS2 (nomodel mode, extsize=200)
# - Generate BED files of accessible regions
# - Create narrowPeak files for downstream analysis
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
FILTER_DIR="$PROJECT_ROOT/results/alignment/filtered"
PEAK_DIR="$PROJECT_ROOT/results/peaks"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$PEAK_DIR" "$LOG_DIR"

# Default parameters (ATAC-seq specific)
GENOME_SIZE="${GENOME_SIZE:-hs}"  # hs = 3.0e9, can also use 3.1e9 or 2.7e9
QVALUE="${QVALUE:-0.01}"
EXT_SIZE="${EXT_SIZE:-200}"  # Nucleosomal spacing for ATAC
SHIFT="${SHIFT:-0}"  # No shift for paired-end BAMPE
NOMODEL=true
FORMAT="${FORMAT:-BAMPE}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

call_peaks_macs2() {
    local input_bam="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Calling peaks with MACS2 for: $sample"
    
    # Check for MACS2
    if ! command -v macs2 &> /dev/null; then
        log "ERROR: MACS2 not found. Install with: pip install macs2"
        return 1
    fi
    
    cd "$output_dir"
    
    # ATAC-seq specific parameters:
    # - nomodel: Don't build shifting model
    # - extsize 200: Fragment size (nucleosomal periodicity)
    # - shift 0: No shift for BAMPE (already have fragment sizes)
    # - format BAMPE: Paired-end mode
    # - g hs: Human genome size (3e9)
    # - q 0.01: q-value threshold
    # - --keep-dup all: Let MACS2 handle duplicates
    
    macs2 callpeak \
        -t "$input_bam" \
        -f "$FORMAT" \
        -g "$GENOME_SIZE" \
        -n "$sample" \
        --nomodel \
        --extsize "$EXT_SIZE" \
        --shift "$SHIFT" \
        -q "$QVALUE" \
        --keep-dup all \
        --outdir "$output_dir" \
        2>&1 | tee -a "$LOG_DIR/macs2_${sample}.log"
    
    cd - > /dev/null
    
    log "Peak calling complete for: $sample"
}

call_peaks_broad() {
    local input_bam="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Calling broad peaks with MACS2 for: $sample"
    
    macs2 callpeak \
        -t "$input_bam" \
        -f "$FORMAT" \
        -g "$GENOME_SIZE" \
        -n "${sample}_broad" \
        --nomodel \
        --extsize "$EXT_SIZE" \
        --shift "$SHIFT" \
        -q 0.05 \
        --broad \
        --broad-cutoff 0.05 \
        --keep-dup all \
        --outdir "$output_dir" \
        2>&1 | tee -a "$LOG_DIR/macs2_${sample}_broad.log"
    
    log "Broad peak calling complete for: $sample"
}

convert_to_bed() {
    local narrowpeak_file="$1"
    local output_bed="$2"
    
    log "Converting narrowPeak to BED format..."
    
    if [[ ! -f "$narrowpeak_file" ]]; then
        log "WARNING: narrowPeak file not found: $narrowpeak_file"
        return 1
    fi
    
    # Convert from narrowPeak (6 mandatory + 10 optional columns) to BED (6 columns)
    # narrowPeak columns: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak
    cut -f1-6 "$narrowpeak_file" > "$output_bed"
    
    log "BED file created: $output_bed"
}

get_summits() {
    local narrowpeak_file="$1"
    local output_bed="$2"
    
    log "Extracting summit positions..."
    
    if [[ ! -f "$narrowpeak_file" ]]; then
        log "WARNING: narrowPeak file not found"
        return 1
    fi
    
    # Get summit positions (peak column = offset from start with highest signal)
    # Add summit offset to start position
    awk 'BEGIN{OFS="\t"} {
        summit = $2 + $10;
        print $1, summit, summit+1, $4, $5, $6
    }' "$narrowpeak_file" > "$output_bed"
    
    log "Summit BED created: $output_bed"
}

merge_peaks() {
    local sample_id="$1"
    local input_dir="$2"
    local output_dir="$3"
    
    log "Merging peak files..."
    
    if ! command -v bedtools &> /dev/null; then
        log "WARNING: bedtools not found. Install with: conda install -c bioconda bedtools"
        return 1
    fi
    
    # Find all narrowPeak files
    local peak_files=$(find "$input_dir" -name "*_peaks.narrowPeak" 2>/dev/null)
    
    if [[ -z "$peak_files" ]]; then
        log "No peak files found to merge"
        return 1
    fi
    
    # Merge peaks from all replicates
    cat $peak_files | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - > "$output_dir/${sample_id}_consensus_peaks.bed"
    
    log "Consensus peaks created: $output_dir/${sample_id}_consensus_peaks.bed"
    echo "$output_dir/${sample_id}_consensus_peaks.bed"
}

peak_annotation() {
    local peak_file="$1"
    local sample="$2"
    local output_dir="$3"
    local annotation_gtf="${ANNOTATION_GTF:-}"
    
    log "Annotating peaks..."
    
    if [[ ! -f "$peak_file" ]]; then
        log "WARNING: Peak file not found"
        return 1
    fi
    
    if ! command -v bedtools &> /dev/null; then
        log "WARNING: bedtools not found, skipping annotation"
        return 1
    fi
    
    # Check for annotation file
    if [[ -z "$annotation_gtf" || ! -f "$annotation_gtf" ]]; then
        log "Annotation GTF not found. Skipping detailed annotation."
        log "Install HOMER for full annotation: conda install -c bioconda homer"
        return 1
    fi
    
    # Simple annotation using bedtools closest
    bedtools closest -a "$peak_file" -b "$annotation_gtf" -D write \
        > "$output_dir/${sample}_annotated_peaks.txt" 2>/dev/null
    
    log "Annotation complete: $output_dir/${sample}_annotated_peaks.txt"
}

generate_peak_report() {
    local sample="$1"
    local peak_dir="$2"
    local output_file="$3"
    
    {
        echo "ATAC-seq Peak Calling Report: $sample"
        echo "====================================="
        echo "Date: $(date)"
        echo ""
        echo "Parameters:"
        echo "  MACS2 version: $(macs2 --version 2>/dev/null || echo 'unknown')"
        echo "  Genome size: $GENOME_SIZE"
        echo "  q-value threshold: $QVALUE"
        echo "  Extension size: $EXT_SIZE"
        echo "  Format: $FORMAT (paired-end)"
        echo ""
        
        local narrowpeak="${peak_dir}/${sample}_peaks.narrowPeak"
        local bed="${peak_dir}/${sample}_peaks.bed"
        local summits="${peak_dir}/${sample}_summits.bed"
        
        echo "Output Files:"
        echo "  narrowPeak: $narrowpeak"
        echo "  BED: $bed"
        echo "  Summits: $summits"
        echo ""
        
        echo "Peak Statistics:"
        echo "----------------"
        if [[ -f "$narrowpeak" ]]; then
            echo "  Total peaks: $(wc -l < "$narrowpeak")"
            echo "  Peaks on chr1: $(grep -c "^chr1	" "$narrowpeak" || echo 0)"
            echo "  Peaks on chrX: $(grep -c "^chrX	" "$narrowpeak" || echo 0)"
            
            echo ""
            echo "SignalValue distribution (top 10):"
            cut -f7 "$narrowpeak" | sort -rn | head -10 | awk '{printf "    %.4f\n", $1}'
            
            echo ""
            echo "q-value distribution (top 10):"
            cut -f9 "$narrowpeak" | sort -rn | head -10 | awk '{printf "    %.2e\n", $1}'
        else
            echo "  Peak files not found"
        fi
        
    } > "$output_file"
    
    log "Peak report saved: $output_file"
}

# Main execution
main() {
    local sample_id="${1:-GM12878_ATAC_Rep1}"
    
    log "=== Step 5: Peak Calling ==="
    
    # Find input BAM from filtering step
    local input_bam="$FILTER_DIR/${sample_id}.filtered.sorted.bam"
    
    if [[ ! -f "$input_bam" ]]; then
        log "ERROR: Input BAM not found: $input_bam"
        log "Please run filtering step first"
        exit 1
    fi
    
    log "Input BAM: $input_bam"
    
    # Step 5a: Call narrow peaks
    call_peaks_macs2 "$input_bam" "$sample_id" "$PEAK_DIR"
    
    # Step 5b: Call broad peaks (optional, for wider regulatory regions)
    call_peaks_broad "$input_bam" "$sample_id" "$PEAK_DIR"
    
    # Step 5c: Convert to BED format
    local narrowpeak="$PEAK_DIR/${sample_id}_peaks.narrowPeak"
    local bed="$PEAK_DIR/${sample_id}_peaks.bed"
    local summits="$PEAK_DIR/${sample_id}_summits.bed"
    
    if [[ -f "$narrowpeak" ]]; then
        convert_to_bed "$narrowpeak" "$bed"
        get_summits "$narrowpeak" "$summits"
    fi
    
    # Step 5d: Merge peaks from replicates (if available)
    merge_peaks "$sample_id" "$PEAK_DIR" "$PEAK_DIR"
    
    # Step 5e: Generate report
    generate_peak_report "$sample_id" "$PEAK_DIR" "$PEAK_DIR/${sample_id}_peak_report.txt"
    
    log "=== Step 5 Complete ==="
    log "Peak files: $PEAK_DIR"
}

main "$@"
