#!/usr/bin/env bash
# ============================================================================
# ATAC-seq Analysis Pipeline - Step 1: Data Acquisition
# ============================================================================
# Downloads public ATAC-seq data from ENCODE or SRA
# Dataset: GM12878 ATAC-seq (ENCODE ENCSR000AAL / GEO GSE65360)
# ============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$PROJECT_ROOT/config"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
LOG_DIR="$PROJECT_ROOT/logs"

# Create directories
mkdir -p "$DATA_DIR" "$LOG_DIR"

# Load configuration
source "$CONFIG_DIR/config.sh" 2>/dev/null || true

# ENCODE GM12878 ATAC-seq download URLs
# ENCODE accession: ENCSR000AAL (GM12878 ATAC-seq, 2 replicates)
ENCODE_URL_1="https://www.encodeproject.org/files/ENCFF234QEM/@@download/ENCFF234QEM.fastq.gz"
ENCODE_URL_2="https://www.encodeproject.org/files/ENCFF235ZKY/@@download/ENCFF235ZKY.fastq.gz"

# Alternative: SRA accessions for GSE65360 (Buenrostro 2018)
SRA_ACCESSIONS=(
    "SRR5626302"  # GM12878 ATAC-seq replicate 1
    "SRR5626303"  # GM12878 ATAC-seq replicate 2
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

download_from_encode() {
    local output_dir="$1"
    local force="$2"
    
    log "Downloading ENCODE GM12878 ATAC-seq data..."
    
    # Check if already downloaded
    if [[ -f "$output_dir/ENCFF234QEM.fastq.gz" && -f "$output_dir/ENCFF235ZKY.fastq.gz" ]]; then
        log "Files already exist. Skipping download (use --force to redownload)"
        [[ "$force" != "true" ]] && return 0
    fi
    
    # Download using curl with retry
    local max_retries=3
    local retry_delay=30
    
    for url in "$ENCODE_URL_1" "$ENCODE_URL_2"; do
        local filename=$(basename "$url")
        local filepath="$output_dir/$filename"
        
        log "Downloading: $filename"
        
        for attempt in $(seq 1 $max_retries); do
            if curl -L -o "$filepath" "$url" 2>/dev/null; then
                log "Successfully downloaded: $filename"
                break
            else
                log "Attempt $attempt/$max_retries failed for $filename"
                [[ $attempt -lt $max_retries ]] && sleep $retry_delay
            fi
        done
        
        # Validate download
        if [[ -f "$filepath" && $(stat -f%z "$filepath" 2>/dev/null || stat -c%s "$filepath" 2>/dev/null) -lt 1000000 ]]; then
            log "ERROR: Downloaded file too small - possible error"
            return 1
        fi
    done
    
    log "ENCODE download complete!"
}

download_from_sra() {
    local accession="$1"
    local output_dir="$2"
    
    log "Downloading from SRA: $accession"
    
    # Check if prefetch and fasterq-dump are available
    if ! command -v prefetch &> /dev/null; then
        log "ERROR: SRA Toolkit (prefetch) not found"
        log "Install with: conda install -c bioconda sra-tools"
        return 1
    fi
    
    # Download SRA file
    prefetch "$accession" -O "$output_dir/sra/"
    
    # Convert to FASTQ
    fasterq-dump "$output_dir/sra/${accession}.sra" -O "$output_dir" --split-files
    
    log "SRA download complete: $accession"
}

create_symlinks() {
    local data_dir="$1"
    
    log "Creating symlinks for pipeline compatibility..."
    
    # Create sample-specific symlinks
    if [[ -f "$data_dir/ENCFF234QEM.fastq.gz" ]]; then
        ln -sf "ENCFF234QEM.fastq.gz" "$data_dir/GM12878_R1.fastq.gz" 2>/dev/null || true
    fi
    
    log "Symlinks created"
}

generate_download_report() {
    local output_file="$1"
    local data_dir="$2"
    
    {
        echo "ATAC-seq Data Download Report"
        echo "============================="
        echo "Date: $(date)"
        echo "Source: ENCODE (ENCSR000AAL) / SRA (SRR5626302)"
        echo ""
        echo "Downloaded Files:"
        ls -lh "$data_dir"/*.fastq.gz 2>/dev/null || echo "No files found"
        echo ""
        echo "Checksums (MD5):"
        find "$data_dir" -name "*.fastq.gz" -exec md5sum {} \; 2>/dev/null || true
    } > "$output_file"
    
    log "Report saved: $output_file"
}

# Main execution
main() {
    local force_download="false"
    local source="encode"
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --force|-f)
                force_download="true"
                shift
                ;;
            --source|-s)
                source="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                exit 1
                ;;
        esac
    done
    
    log "=== Step 1: Data Acquisition ==="
    log "Source: ${source}"
    
    case "$source" in
        encode)
            download_from_encode "$DATA_DIR" "$force_download"
            ;;
        sra)
            for acc in "${SRA_ACCESSIONS[@]}"; do
                download_from_sra "$acc" "$DATA_DIR"
            done
            ;;
        *)
            log "Unknown source: $source"
            exit 1
            ;;
    esac
    
    create_symlinks "$DATA_DIR"
    generate_download_report "$RESULTS_DIR/download_report.txt" "$DATA_DIR"
    
    log "=== Step 1 Complete ==="
    log "Raw data location: $DATA_DIR"
}

main "$@"
