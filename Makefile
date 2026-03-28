# ATAC-seq Pipeline Makefile
# Simplified workflow execution

.PHONY: help install download qc trim align filter peaks visualize all clean

# Configuration
SAMPLE ?= GM12878_ATAC_Rep1
THREADS ?= 8
CONFIG = config/config.yaml

# Directories
DATA = data
RESULTS = results

help:
	@echo "ATAC-seq Analysis Pipeline"
	@echo "========================="
	@echo ""
	@echo "Available targets:"
	@echo "  install    - Create conda environments"
	@echo "  download   - Download ENCODE GM12878 ATAC-seq data"
	@echo "  qc         - Run FastQC on raw reads"
	@echo "  trim       - Trim adapters and low-quality bases"
	@echo "  align      - Align reads with Bowtie2"
	@echo "  filter     - Remove duplicates, filter by MAPQ"
	@echo "  peaks      - Call peaks with MACS2"
	@echo "  visualize  - Generate R visualizations"
	@echo "  all        - Run complete pipeline"
	@echo "  clean      - Remove all results"
	@echo ""
	@echo "Usage:"
	@echo "  make download SAMPLE=SRR5626302"
	@echo "  make align THREADS=16"
	@echo "  make all"

install:
	@echo "Creating conda environments..."
	conda env create -f envs/atacseq.yaml
	conda env create -f envs/r-vis.yaml
	@echo "Installation complete!"

download:
	@echo "Downloading ATAC-seq data..."
	bash scripts/01_download.sh

qc:
	@echo "Running FastQC..."
	bash scripts/02_qc_trim.sh --skip-trim

trim:
	@echo "Trimming adapters..."
	bash scripts/02_qc_trim.sh

align:
	@echo "Aligning reads..."
	bash scripts/03_align.sh

filter:
	@echo "Filtering BAM files..."
	bash scripts/04_filter.sh

peaks:
	@echo "Calling peaks..."
	bash scripts/05_peak_calling.sh

footprints:
	@echo "Running motif analysis..."
	bash scripts/06_footprinting.sh

visualize:
	@echo "Generating visualizations..."
	Rscript scripts/07_visualization.R

multiqc:
	@echo "Generating MultiQC report..."
	multiqc results/fastqc -o results -n multiqc_report

all: download trim align filter peaks visualize multiqc
	@echo ""
	@echo "Pipeline complete! Check results/ directory."

clean:
	@echo "Cleaning results..."
	rm -rf $(DATA)/*.fastq.gz
	rm -rf results/fastqc/*.html results/fastqc/*.zip
	rm -rf results/trimming/*.fastq.gz results/trimming/*.txt
	rm -rf results/alignment/*.bam results/alignment/*.bai results/alignment/*.sam
	rm -rf results/alignment/filtered/*.bam results/alignment/filtered/*.bai
	rm -rf results/peaks/*
	rm -rf results/footprints/*
	rm -rf results/plots/*.png
	rm -f results/multiqc_report.html
	@echo "Clean complete!"

# Snakemake integration
snakemake:
	snakemake -p --use-conda --cores $(THREADS) -s SnakeMake/Snakefile

dry-run:
	snakemake -n -p -s SnakeMake/Snakefile

# Development
test:
	@echo "Testing pipeline setup..."
	@bash scripts/01_download.sh --help || true
	@Rscript --version || true
