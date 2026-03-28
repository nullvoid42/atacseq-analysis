# ATAC-seq Chromatin Accessibility Analysis Pipeline

A complete, reproducible ATAC-seq analysis pipeline for characterizing chromatin accessibility, 
reproducing the classic Buenrostro 2018 (GEO: GSE65360) case study using GM12878 lymphoblastoid cells.

## 🧬 Overview

This pipeline implements the full ATAC-seq analysis workflow from raw sequencing reads to 
transcription factor footprinting and publication-quality visualizations.

**Primary Dataset:** ENCODE GM12878 ATAC-seq (ENCFF234QEM)  
**Reference:** Buenrostro JD, et al. (2018) *Nature Methods* | GEO: GSE65360

## 📁 Project Structure

```
atacseq-analysis/
├── config/
│   ├── config.yaml           # Main configuration
│   └── samples.tsv           # Sample manifest
├── scripts/
│   ├── 01_download.sh        # Data acquisition
│   ├── 02_qc_trim.sh         # FastQC + adapter trimming
│   ├── 03_align.sh           # Bowtie2 alignment
│   ├── 04_filter.sh          # Deduplication & filtering
│   ├── 05_peak_calling.sh    # MACS2 peak calling
│   ├── 06_footprinting.sh    # TF motif annotation
│   └── 07_visualization.R     # Comprehensive visualizations
├── SnakeMake/
│   └── Snakefile             # Snakemake workflow
├── results/
│   ├── fastqc/               # FastQC reports
│   ├── trimming/             # Trimmed reads
│   ├── alignment/            # BAM files
│   ├── peaks/                # MACS2 outputs
│   ├── footprints/           # Motif analysis
│   └── plots/                # Visualizations
├── envs/                     # Conda environments
│   ├── atacseq.yaml
│   └── r-vis.yaml
├── CITATION.cff
└── README.md
```

## 🔧 Installation

### Prerequisites

```bash
# Install Conda/Mamba
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

# Clone this repository
git clone https://github.com/nullvoid42/atacseq-analysis.git
cd atacseq-analysis

# Create environments
conda env create -f envs/atacseq.yaml
conda env create -f envs/r-vis.yaml
```

### Required Tools

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | ≥0.11.9 | Read quality assessment |
| Cutadapt | ≥4.0 | Adapter trimming |
| Bowtie2 | ≥2.5.0 | Read alignment |
| SAMtools | ≥1.17 | BAM processing |
| Picard | ≥3.0.0 | Duplicate removal |
| MACS2 | ≥2.7.0 | Peak calling |
| bedtools | ≥2.30.0 | BED operations |
| deepTools | ≥3.5.0 | BigWig & heatmaps |
| HOMER | ≥4.11 | Motif annotation |
| R | ≥4.2.0 | Visualization |

## 🚀 Quick Start

### 1. Configure your analysis

Edit `config/config.yaml`:

```yaml
project: GM12878_ATACseq
genome: hg38
reference: /path/to/bowtie2_index/hg38
samples:
  - id: GM12878_ATAC
    fastq_1: /data/GM12878_R1.fastq.gz
    fastq_2: /data/GM12878_R2.fastq.gz
```

### 2. Download test data (optional)

```bash
# Download from ENCODE
bash scripts/01_download.sh

# Or from SRA
prefetch SRR12345678
fasterq-dump SRR12345678 --split-files
```

### 3. Run the pipeline

**Using Snakemake (recommended):**
```bash
conda activate atacseq
snakemake -p --use-conda --cores 8

# Dry run first
snakemake -n -p
```

**Using individual scripts:**
```bash
bash scripts/01_download.sh
bash scripts/02_qc_trim.sh
bash scripts/03_align.sh
bash scripts/04_filter.sh
bash scripts/05_peak_calling.sh
bash scripts/06_footprinting.sh
Rscript scripts/07_visualization.R
```

## 📊 Pipeline Overview

### Step 1: Data Acquisition
- Download ENCODE ATAC-seq GM12878 data (ENCFF234QEM)
- SRA toolkit for NCBI data retrieval
- Checksum validation

### Step 2: Quality Control & Trimming
- **FastQC**: Per-base quality, adapter content, read length distribution
- **Cutadapt**: Remove Nextera XT adapters (CTGTCTCTTATACACATCT)
- **chrM removal**: Filter mitochondrial reads

### Step 3: Alignment
- **Bowtie2**: very-sensitive mode, paired-end
- **SAMtools**: Convert, sort, index (BAM)
- Reference: GRCh38/hg38

### Step 4: Filtering
- **Picard MarkDuplicates**: Remove PCR duplicates
- **SAMtools filter**: MAPQ ≥ 30, proper pairs
- **chrM exclusion**: Remove mitochondrial alignments

### Step 5: Peak Calling
- **MACS2 callpeak**: nomodel, extsize=200, shift=0
- Generates: peaks.bed, summits.bed, narrowPeak
- Automatic q-value thresholding

### Step 6: Footprinting (Optional)
- **HOMER findMotifsGenome**: De novo motif discovery
- ** annotatePeaks.pl**: Known motif annotation
- **PWMScan**: TF binding site identification

### Step 7: Visualization
- **BigWig tracks**: Normalized coverage (RPKM)
- **TSS heatmaps**: Accessibility around transcription start sites
- **Peak annotations**: Genomic distribution (promoter, intron, intergenic)
- **Motif logos**: Top enriched TFs
- **FRiP score**: Fraction of reads in peaks

## 📈 Expected Results

### Quality Metrics

| Metric | Expected Value |
|--------|---------------|
| Raw reads | 50-100M (ENCODE GM12878) |
| Fragment size | ~200bp median (nucleosomal pattern) |
| Unique mapping rate | >70% |
| FRiP score | >0.30 |
| Number of peaks | 30,000-100,000 |

### Key Findings (Buenrostro 2018 GM12878)

1. **Open chromatin regions**: ~50,000-80,000 accessible sites
2. **Peak distribution**: Prominent enrichment at promoters and enhancers
3. **Top TFs**: CTCF, NFKB, RELA, ETS family motifs
4. **Nucleosomal pattern**: Clear 1-3 nucleosome banding in fragment size distribution

## 🔬 Methods Summary

### Alignment Parameters
```bash
bowtie2 --very-sensitive \
        -X 2000 \
        --fr \
        -x hg38 \
        -1 R1.fastq.gz \
        -2 R2.fastq.gz \
| samtools view -bS \
| samtools sort -o aligned.bam
```

### Peak Calling Parameters
```bash
macs2 callpeak \
    -t aligned.filtered.bam \
    -f BAMPE \
    -g hs \
    -n GM12878_ATAC \
    --nomodel \
    --extsize 200 \
    -q 0.01
```

### Normalization
- **BigWig**: RPKM normalization (1x scaling)
- **Heatmaps**: CPM normalization, Z-score row-wise
- **Peak calling**: Dynamic background estimation

## 📚 References

- Buenrostro JD, et al. (2018) Single chromatin accessibility reveals 
  multi-kinase and transcriptional regulation of hematopoiesis. 
  *Nature Methods* 15: 1006-1012. GEO: GSE65360

- ENCODE Project Consortium (2012) An integrated encyclopedia of DNA 
  elements in the human genome. *Nature* 489: 57-74.

- Buenrostro JD, et al. (2013) Transposition of native chromatin for 
  fast and sensitive epigenomic profiling. *Nature Methods* 10: 1213-1218.

## 📝 License

MIT License - see LICENSE file.

## 🙋 Support

Open an issue on GitHub for questions or problems.
