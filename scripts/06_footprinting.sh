#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 6: Transcription Factor Footprinting
# ============================================================================
# - Motif annotation with HOMER
# - De novo motif discovery
# - Known motif annotation
# - TF binding site identification
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
PEAK_DIR="$PROJECT_ROOT/results/peaks"
FOOTPRINT_DIR="$PROJECT_ROOT/results/footprints"
LOG_DIR="$PROJECT_ROOT/logs"
GENOME="${GENOME:-hg38}"

mkdir -p "$FOOTPRINT_DIR" "$LOG_DIR"

# HOMER parameters
MOTIF_LENGTHS="${MOTIF_LENGTHS:-6,8,10,12}"
ORGANISM="${ORGANISM:-hg38}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

find_motifs_genome() {
    local peak_file="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Finding motifs with HOMER for: $sample"
    
    # Check for HOMER
    if ! command -v findMotifsGenome.pl &> /dev/null; then
        log "WARNING: HOMER not found. Install with: conda install -c bioconda homer"
        log "Falling back to basic motif analysis..."
        basic_motif_analysis "$peak_file" "$sample" "$output_dir"
        return 1
    fi
    
    local motif_output="$output_dir/${sample}_motifs"
    mkdir -p "$motif_output"
    
    # HOMER findMotifsGenome for de novo motif discovery
    findMotifsGenome.pl \
        "$peak_file" \
        "$ORGANISM" \
        "$motif_output" \
        -length "$MOTIF_LENGTHS" \
        -size given \
        -mis 2 \
        2>&1 | tee -a "$LOG_DIR/homer_${sample}.log"
    
    log "De novo motif analysis complete: $motif_output"
}

annotate_peaks() {
    local peak_file="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Annotating peaks with HOMER..."
    
    if ! command -v annotatePeaks.pl &> /dev/null; then
        log "WARNING: HOMER annotatePeaks.pl not found"
        return 1
    fi
    
    local annotated="$output_dir/${sample}_annotated.tsv"
    
    annotatePeaks.pl \
        "$peak_file" \
        "$ORGANISM" \
        > "$annotated" \
        2>&1 | tee -a "$LOG_DIR/homer_annotate_${sample}.log"
    
    log "Peak annotation complete: $annotated"
}

basic_motif_analysis() {
    local peak_file="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Running basic motif analysis (without HOMER)..."
    
    # Create a simple motif report without HOMER
    {
        echo "Motif Analysis Report: $sample"
        echo "============================="
        echo "Date: $(date)"
        echo "Note: HOMER not installed - basic analysis only"
        echo ""
        echo "Top known TF motifs expected in ATAC-seq data:"
        echo "  1. CTCF - CCCTC-binding factor (insulator)"
        echo "  2. NF-kB (p65/RELA) - Inflammatory response"
        echo "  3. AP-1 (FOS/JUN) - Cell proliferation"
        echo "  4. TEAD - Hippo signaling"
        echo "  5. ETS family (ELK1, GABPA) - Growth factor response"
        echo "  6. MYC - Oncogene"
        echo "  7. MAX - Cell growth"
        echo "  8. RFX - Immune response"
        echo ""
        echo "For full motif analysis, install HOMER:"
        echo "  conda install -c bioconda homer"
        echo ""
        echo "Peak file: $peak_file"
        echo "Peak count: $(wc -l < "$peak_file" 2>/dev/null || echo 'N/A')"
        
    } > "$output_dir/${sample}_motif_report.txt"
}

scan_motifs() {
    local summit_bed="$1"
    local sample="$2"
    local output_dir="$3"
    
    log "Scanning for known TF motifs..."
    
    # Check for HOMER
    if ! command -v findMotifsGenome.pl &> /dev/null; then
        log "HOMER not found - skipping motif scan"
        return 1
    fi
    
    # Known motif scanning
    findMotifsGenome.pl \
        "$summit_bed" \
        "$ORGANISM" \
        "$output_dir/${sample}_known_motifs" \
        -length "$MOTIF_LENGTHS" \
        -size 200 \
        -海量搜索known \
        2>&1 | tee -a "$LOG_DIR/homer_known_${sample}.log"
    
    log "Known motif scanning complete"
}

generate_footprint_report() {
    local sample="$1"
    local footprint_dir="$2"
    local output_file="$3"
    
    {
        echo "ATAC-seq Footprinting Report: $sample"
        echo "======================================"
        echo "Date: $(date)"
        echo ""
        echo "Analysis Method: HOMER findMotifsGenome.pl"
        echo "Organism: $ORGANISM"
        echo "Motif lengths: $MOTIF_LENGTHS"
        echo ""
        
        local motif_dir="$footprint_dir/${sample}_motifs"
        if [[ -d "$motif_dir" ]]; then
            echo "De novo Motif Results:"
            echo "----------------------"
            if [[ -f "$motif_dir"/knownResults.html ]]; then
                echo "Known motif results: $motif_dir/knownResults.html"
            fi
            if [[ -f "$motif_dir"/homerResults.html ]]; then
                echo "De novo results: $motif_dir/homerResults.html"
            fi
        fi
        
        echo ""
        echo "Expected TF motifs in GM12878 (Lymphoblastoid cells):"
        echo "  - CTCF (insulator protein)"
        echo "  - NF-kB/RELA (EBV latent infection marker)"
        echo "  - IRF4 (immune response)"
        echo "  - EBF1 (B cell development)"
        echo "  - PAX5 (B cell specification)"
        
    } > "$output_file"
    
    log "Footprint report saved: $output_file"
}

# Main execution
main() {
    local sample_id="${1:-GM12878_ATAC_Rep1}"
    
    log "=== Step 6: TF Footprinting ==="
    
    # Find input peak files
    local peak_bed="$PEAK_DIR/${sample_id}_peaks.bed"
    local summit_bed="$PEAK_DIR/${sample_id}_summits.bed"
    
    if [[ ! -f "$peak_bed" ]]; then
        log "ERROR: Peak BED not found: $peak_bed"
        log "Please run peak calling step first"
        exit 1
    fi
    
    log "Input peaks: $peak_bed"
    log "Input summits: $summit_bed"
    
    # Step 6a: De novo motif discovery
    find_motifs_genome "$peak_bed" "$sample_id" "$FOOTPRINT_DIR"
    
    # Step 6b: Annotate peaks
    annotate_peaks "$peak_bed" "$sample_id" "$FOOTPRINT_DIR"
    
    # Step 6c: Scan for known motifs
    if [[ -f "$summit_bed" ]]; then
        scan_motifs "$summit_bed" "$sample_id" "$FOOTPRINT_DIR"
    fi
    
    # Step 6d: Generate report
    generate_footprint_report "$sample_id" "$FOOTPRINT_DIR" \
        "$FOOTPRINT_DIR/${sample_id}_footprint_report.txt"
    
    log "=== Step 6 Complete ==="
    log "Footprint results: $FOOTPRINT_DIR"
}

main "$@"
