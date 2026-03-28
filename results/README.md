# ATAC-seq Analysis Results

This directory will contain the analysis outputs after running the pipeline.

## Directory Structure

```
results/
├── raw_data/           # Downloaded FASTQ files
├── fastqc/             # FastQC quality reports
│   ├── raw/           # Raw read QC
│   └── trimmed/       # Trimmed read QC
├── trimming/          # Adapter-trimmed FASTQ files
├── alignment/         # Bowtie2 alignment outputs
│   └── filtered/      # Filtered BAM files (MAPQ >= 30, no chrM)
├── peaks/             # MACS2 peak calling results
│   ├── *_peaks.narrowPeak
│   ├── *_peaks.bed
│   └── *_summits.bed
├── footprints/        # HOMER motif analysis
├── plots/             # R visualization outputs
│   ├── 01_fragment_size_distribution.png
│   ├── 02_peak_annotation_pie.png
│   ├── 03_tss_heatmap.png
│   └── ...
└── multiqc_report.html
```

## Expected Outputs

### Quality Control
- `fastqc/*_fastqc.html` - Per-sample quality reports
- `multiqc_report.html` - Aggregated QC report

### Alignment
- `alignment/*.sorted.bam` - Sorted BAM files
- `alignment/*.sorted.bam.bai` - BAM indices

### Filtering
- `alignment/filtered/*.filtered.sorted.bam` - High-quality filtered BAM
- `alignment/filtered/*.dup_metrics.txt` - Duplicate statistics
- `alignment/filtered/*.fragment_sizes.txt` - Fragment size distribution

### Peak Calling
- `peaks/*_peaks.narrowPeak` - NarrowPeak format (MACS2)
- `peaks/*_peaks.bed` - BED format for genome browsers
- `peaks/*_summits.bed` - Peak summit positions
- `peaks/*_peak_report.txt` - Peak calling statistics

### Footprinting
- `footprints/*_motifs/` - HOMER motif discovery results
- `footprints/*_annotated.tsv` - Peak annotations

### Visualization
- `plots/01_fragment_size_distribution.png` - Nucleosomal pattern
- `plots/02_peak_annotation_pie.png` - Genomic annotation pie chart
- `plots/03_tss_heatmap.png` - TSS accessibility heatmap
- `plots/04_peak_width_distribution.png` - Peak width density
- `plots/05_genomic_distribution.png` - Genomic region barplot
- `plots/06_peak_overlap.png` - Replicate overlap heatmap
- `plots/07_quality_dashboard.png` - Quality metrics summary
- `plots/08_summary_statistics.png` - Overall analysis summary

## Key Metrics to Check

| Metric | Expected Value | Interpretation |
|--------|---------------|----------------|
| Unique mapping rate | >70% | Good library complexity |
| FRiP score | >0.30 | Significant open chromatin |
| Number of peaks | 30,000-100,000 | Typical for GM12878 |
| Median fragment size | ~200bp | Mononucleosomal peak |
| Promoter peaks | 30-40% | Active gene regulation |
