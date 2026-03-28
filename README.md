# ATAC-seq 染色质可及性分析流程

一套完整、可复现的 ATAC-seq 分析流程，用于表征染色质可及性，基于经典 Buetrostro 2018（GEO：GSE65360）案例研究，使用 GM12878 淋巴母细胞系。

## 🧬 概述

本流程实现了从原始测序读数到转录因子足迹分析和出版级可视化的完整 ATAC-seq 分析工作流程。

**主要数据集：** ENCODE GM12878 ATAC-seq（ENCFF234QEM）
**参考文献：** Buenrostro JD, et al. (2018) *Nature Methods* | GEO: GSE65360

## 📁 项目结构

```
atacseq-analysis/
├── config/
│   ├── config.yaml           # 主配置文件
│   └── samples.tsv           # 样本清单
├── scripts/
│   ├── 01_download.sh        # 数据获取
│   ├── 02_qc_trim.sh         # FastQC 质控 + 接头剪切
│   ├── 03_align.sh           # Bowtie2 比对
│   ├── 04_filter.sh          # 去重与过滤
│   ├── 05_peak_calling.sh    # MACS2 峰值检测
│   ├── 06_footprinting.sh    # TF 基序注释
│   └── 07_visualization.R     # 综合可视化
├── SnakeMake/
│   └── Snakefile             # Snakemake 工作流
├── results/
│   ├── fastqc/               # FastQC 报告
│   ├── trimming/             # 剪切后reads
│   ├── alignment/            # BAM 文件
│   ├── peaks/                # MACS2 输出
│   ├── footprints/           # 基序分析
│   └── plots/                # 可视化结果
├── envs/                     # Conda 环境
│   ├── atacseq.yaml
│   └── r-vis.yaml
├── CITATION.cff
└── README.md
```

## 🔧 安装

### 前置条件

```bash
# 安装 Conda/Mamba
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

# 克隆本仓库
git clone https://github.com/nullvoid42/atacseq-analysis.git
cd atacseq-analysis

# 创建环境
conda env create -f envs/atacseq.yaml
conda env create -f envs/r-vis.yaml
```

### 所需工具

| 工具 | 版本 | 用途 |
|------|------|------|
| FastQC | ≥0.11.9 | 读数质量评估 |
| Cutadapt | ≥4.0 | 接头剪切 |
| Bowtie2 | ≥2.5.0 | 序列比对 |
| SAMtools | ≥1.17 | BAM 处理 |
| Picard | ≥3.0.0 | 去重 |
| MACS2 | ≥2.7.0 | 峰值检测 |
| bedtools | ≥2.30.0 | BED 操作 |
| deepTools | ≥3.5.0 | BigWig 和热图 |
| HOMER | ≥4.11 | 基序注释 |
| R | ≥4.2.0 | 可视化 |

## 🚀 快速开始

### 1. 配置分析参数

编辑 `config/config.yaml`：

```yaml
project: GM12878_ATACseq
genome: hg38
reference: /path/to/bowtie2_index/hg38
samples:
  - id: GM12878_ATAC
    fastq_1: /data/GM12878_R1.fastq.gz
    fastq_2: /data/GM12878_R2.fastq.gz
```

### 2. 下载测试数据（可选）

```bash
# 从 ENCODE 下载
bash scripts/01_download.sh

# 或从 SRA 下载
prefetch SRR12345678
fasterq-dump SRR12345678 --split-files
```

### 3. 运行流程

**使用 Snakemake（推荐）：**
```bash
conda activate atacseq
snakemake -p --use-conda --cores 8

# 先预览
snakemake -n -p
```

**使用独立脚本：**
```bash
bash scripts/01_download.sh
bash scripts/02_qc_trim.sh
bash scripts/03_align.sh
bash scripts/04_filter.sh
bash scripts/05_peak_calling.sh
bash scripts/06_footprinting.sh
Rscript scripts/07_visualization.R
```

## 📊 流程概览

### 第一步：数据获取
- 下载 ENCODE ATAC-seq GM12878 数据（ENCFF234QEM）
- 使用 SRA toolkit 从 NCBI 获取数据
- 校验和验证

### 第二步：质控与剪切
- **FastQC**：逐碱基质量、接头含量、读长分布
- **Cutadapt**：去除 Nextera XT 接头（CTGTCTCTTATACACATCT）
- **chrM 去除**：过滤线粒体 reads

### 第三步：比对
- **Bowtie2**：超敏感模式，双端测序
- **SAMtools**：转换、排序、建立索引（BAM）
- 参考基因组：GRCh38/hg38

### 第四步：过滤
- **Picard MarkDuplicates**：去除 PCR 重复
- **SAMtools filter**：MAPQ ≥ 30，正确配对
- **chrM 排除**：去除线粒体比对结果

### 第五步：峰值检测
- **MACS2 callpeak**：nomodel，extsize=200，shift=0
- 生成：peaks.bed、summits.bed、narrowPeak
- 自动 q 值阈值化

### 第六步：足迹分析（可选）
- **HOMER findMotifsGenome**：从头发现基序
- **annotatePeaks.pl**：已知基序注释
- **PWMScan**：转录因子结合位点鉴定

### 第七步：可视化
- **BigWig 轨迹**：归一化覆盖度（RPKM）
- **TSS 热图**：转录起始位点周围的可及性
- **峰值注释**：基因组分布（启动子、内含子、基因间区）
- **基序 logo**：富集 top 转录因子
- **FRiP 评分**：落在峰值内的 reads 比例

## 📈 预期结果

### 质控指标

| 指标 | 预期值 |
|------|--------|
| 原始 reads 数 | 50-100M（ENCODE GM12878） |
| 片段长度 | 中位数约 200bp（核小体模式） |
| 唯一比对率 | >70% |
| FRiP 评分 | >0.30 |
| 峰值数量 | 30,000-100,000 |

### 主要发现（Buenrostro 2018 GM12878）

1. **开放染色质区域**：约 50,000-80,000 个可及位点
2. **峰值分布**：启动子和增强子区域显著富集
3. **主要转录因子**：CTCF、NFKB、RELA、ETS 家族基序
4. **核小体模式**：片段长度分布中清晰的 1-3 核小体条带

## 🔬 方法摘要

### 比对参数
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

### 峰值检测参数
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

### 归一化方法
- **BigWig**：RPKM 归一化（1x 缩放）
- **热图**：CPM 归一化，行向 Z-score
- **峰值检测**：动态背景估计

## 📚 参考文献

- Buenrostro JD, et al. (2018) Single chromatin accessibility reveals
  multi-kinase and transcriptional regulation of hematopoiesis.
  *Nature Methods* 15: 1006-1012. GEO: GSE65360

- ENCODE Project Consortium (2012) An integrated encyclopedia of DNA
  elements in the human genome. *Nature* 489: 57-74.

- Buenrostro JD, et al. (2013) Transposition of native chromatin for
  fast and sensitive epigenomic profiling. *Nature Methods* 10: 1213-1218.

## 📝 许可证

MIT License - 详见 LICENSE 文件。

## 🙋 支持

如有疑问或问题，请在 GitHub 上提交 issue。

