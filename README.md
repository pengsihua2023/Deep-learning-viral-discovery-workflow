# Machine Learning-Enhanced Metagenomic Viral Discovery Workflow

**Version 5.2.1** | ML-Powered Novel Virus Discovery | Multi-Tool Validation

A comprehensive Nextflow workflow for discovering and validating viral sequences from metagenomic data using machine learning and Pfam protein validation.

**Key Features**:
- **Machine Learning**: DeepVirFinder for novel virus discovery
- **Pfam Validation**: viralFlye for protein domain verification
- **Multi-Tool Consensus**: Cross-validation for high-confidence identification
- **Viral Abundance Analysis**: RPM and RPKM calculation for all identified viruses

**Supports two modes**:
- **Short-read mode**: Dual assemblers (MEGAHIT + SPAdes) + ML tools (VirSorter2 + DeepVirFinder)
- **Long-read mode**: Three-tool parallel analysis (VirSorter2 + DeepVirFinder + viralFlye) 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.7-blue.svg)](https://www.python.org/)

---

## Table of Contents

- [Core Features](#-core-features)
- [Mode Comparison](#-mode-comparison-summary) 
- [Quick Start](#-quick-start)
- [Three-Tool Design Philosophy](#-three-tool-design-philosophy-long-read-mode)
- [Machine Learning for Novel Virus Discovery](#-machine-learning-for-novel-virus-discovery) 
- [Installation Requirements](#-installation-requirements)
- [Input Data Format](#-input-data-format)
- [Usage](#-usage)
- [Parameter Configuration](#️-parameter-configuration)
- [Output Structure](#-output-structure)
- [Results Interpretation](#-results-interpretation)
- [Performance Optimization](#-performance-and-resources)
- [FAQ](#-frequently-asked-questions)
- [Citation](#-citation)
- **[📘 Mode Selection Guide](./Mode_Selection_Guide.md)**  

                     
---

## Core Features

### 🔬 Dual-Mode Viral Identification Design

#### Long-read Mode: Three-Tool Parallel Analysis 

```
metaFlye assembly
     ↓
 All contigs
     ↓
  ┌──┼──┐
  ↓  ↓  ↓
 VS2 DVF viralFlye
  |  |  |
  └──┼──┘
     ↓
Three-tool comparison
  - 3-tool consensus ⭐⭐⭐ (Highest confidence)
  - 2-tool consensus ⭐⭐ (Medium confidence)
  - Single-tool ⭐ (Exploratory)
```

**Features**:
- ✅ Three tools run **independently in parallel**
- ✅ viralFlye provides **Pfam protein validation**
- ✅ Confidence stratification with clear results

#### Short-read Mode: Dual Assemblers + Dual Tools

```
Illumina Reads
     ↓
  ┌──┴──┐
  ↓     ↓
MEGAHIT SPAdes
  ↓     ↓
  ├─VS2 ├─VS2
  └─DVF └─DVF
  ↓     ↓
  └──┬──┘
     ↓
Assembler comparison
```

**Features**:
- ✅ Dual assemblers in **parallel** (MEGAHIT + SPAdes)
- ✅ Dual-tool verification (VirSorter2 + DeepVirFinder)
- ✅ Cross-assembler comparison, identify consensus viruses

### 🛠️ Viral Identification Tools (ML-Enhanced)

| Tool | Method | Characteristics | Mode Support | Novel Virus Discovery |
|------|--------|----------------|-------------|----------------------|
| **VirSorter2** | Viral DB + Machine Learning | Balanced | Short+Long | Known + Novel  |
| **DeepVirFinder** | Deep Neural Network | High sensitivity | Short+Long | **Excellent for Novel** |
| **viralFlye** | Pfam protein validation | High specificity | **Long-read only** | Function-based Novel Discovery |

**Machine Learning Advantage** :
- **DeepVirFinder**: Learns sequence patterns, discovers viruses **without relying on sequence similarity**
- **VirSorter2**: Combines ML with viral features, identifies novel variants
- **viralFlye**: Function-based approach, discovers viruses **by protein domains** (not sequence)

**Note**: viralFlye **only supports long-read mode** (requires complete metaFlye output)

### Multi-Assembler Support

- **Short-read**: MEGAHIT + metaSPAdes dual assemblers in parallel
- **Long-read**: metaFlye (Nanopore/PacBio)

### Other Features

- ✅ Complete QC pipeline (fastp)
- ✅ Save complete metaFlye output (assembly_info.txt, assembly_graph.gfa)
- ✅ SLURM cluster optimization
- ✅ Apptainer/Singularity container support
- ✅ Optimized parameter configuration (balance sensitivity and specificity)

---

## Quick Start

### 30-Second Quick Launch

```bash
# 1. Prepare sample sheet
echo "sample,fastq_long" > samplesheet_long.csv
echo "sample1,/path/to/reads.fastq.gz" >> samplesheet_long.csv

# 2. Run workflow
sbatch run_metagenome_assembly_classification_longread.sh

# 3. View results (after completion)
cat results_long/three_tools_comparison/*_comparison.txt
```

**That's it!** 

---

## Mode Comparison Summary

### Long-read vs Short-read

| Feature | Short-read Mode | Long-read Mode |
|---------|----------------|----------------|
| **Sequencing Platform** | Illumina | Nanopore / PacBio |
| **Assemblers** | MEGAHIT + metaSPAdes (dual parallel) | metaFlye |
| **Viral Identification** | VirSorter2 + DeepVirFinder (dual) | VirSorter2 + DeepVirFinder + viralFlye (three)  |
| **Specialty** | Dual-assembler cross-validation | **Pfam protein validation** (viralFlye)  |
| **Comparison** | Assembler comparison (MEGAHIT vs SPAdes) | **Three-tool comparison** (VS2 vs DVF vs viralFlye)  |
| **Confidence Tiers** | Dual-tool consensus | **1-3 tool consensus** (more detailed)  |
| **Virus Count** | Medium | More (three tools cover wider range) |
| **viralFlye Support** | ❌ Not supported | ✅ Supported ⭐ |

**Summary**:
- **Short-read mode**: Dual-assembler strategy, focus on assembler consistency
- **Long-read mode**: Three-tool strategy, focus on multi-method validation + Pfam protein validation 

**Detailed comparison**: See [`Mode Selection Guide`](./Mode_Selection_Guide.md) ⭐

---

## Three-Tool Design Philosophy (Long-read Mode)

### Why "Parallel" Instead of "Sequential"?

#### ❌ Old Design (Sequential Validation)

```
metaFlye → VS2+DVF → viralFlye → VS2+DVF validation again
                                    ↑
                                Redundant! Waste of resources
```

#### ✅ New Design (Parallel Independent) 

```
metaFlye → [VS2 ∥ DVF ∥ viralFlye] → Comprehensive comparison
            Run independently, validate each other
```

**Advantages**:
- Scientifically sound (three independent lines of evidence)
- Avoid redundancy (save ~40 hours of computation)
- Clear results (stratified by consensus level)

### Three Tools' Positioning (ML-Enhanced Novel Virus Discovery) 

**VirSorter2** (Hybrid ML):
- **Method**: Viral feature database + **Machine Learning classifiers**
- **Novel Virus Discovery**: Can identify novel variants of known virus families
- **Characteristics**: Balance sensitivity and specificity
- **Output**: Viral classification (dsDNA phage, ssDNA, RNA virus, etc.)

**DeepVirFinder** (Deep Neural Network) :
- **Method**: **Deep Learning** model trained on viral vs non-viral sequences
- **Novel Virus Discovery**: **Excellent** - learns sequence patterns, not similarity
- **Key Advantage**: Discovers viruses **without requiring homology to known viruses**
- **Characteristics**: High sensitivity, pattern-based identification
- **Output**: Viral probability score and p-value

**viralFlye** (Function-Based) :
- **Method**: **Pfam protein domain** validation (function-based)
- **Novel Virus Discovery**: Discovers viruses **by protein function**, not sequence
- **Key Advantage**: Can identify completely novel viruses if they have viral protein domains
- **Characteristics**: Strict validation, high specificity, low false positives
- **Output**: Pfam-validated viruses (Virus/Chromosome/Plasmid classification)

**Why Three Tools for Novel Virus Discovery?** 

```
DeepVirFinder (ML) → Discovers by sequence patterns 
     +
VirSorter2 (Hybrid) → Discovers by features + ML 
     +
viralFlye (Function) → Validates by protein domains 
     ↓
Comprehensive Novel Virus Discovery 
```

**Important**: Pfam is NOT a virus-specific database, but a **universal protein family database** (contains proteins from all life forms). viralFlye distinguishes viruses by analyzing domain combination patterns, enabling **function-based novel virus discovery**.

---

## Machine Learning for Novel Virus Discovery

### Why Machine Learning is Essential for Discovering New Viruses

**Traditional Methods (Limitations)** ❌:
- **BLAST-based**: Requires sequence similarity to known viruses
- **Homology search**: Misses completely novel viruses
- **Marker genes**: Limited to specific virus families

**Machine Learning Approach (This Workflow)** ✅:

#### 1. DeepVirFinder - Deep Neural Network 

```
Trained on thousands of viral and non-viral sequences
    ↓
Learns "what makes a sequence viral"
    ↓
Pattern recognition (NOT sequence similarity)
    ↓
Discovers completely novel viruses 
```

**Key Advantages**:
- ✅ **No homology required** - discovers viruses with 0% similarity to known viruses
- ✅ **Pattern-based** - recognizes viral sequence composition patterns
- ✅ **High sensitivity** - casts a wide net for candidates
- ✅ **P-value confidence** - quantifies discovery confidence

**Example**: Can discover a virus from a completely unexplored virus family!

#### 2. VirSorter2 - Hybrid ML Approach 

```
Viral features (hallmark genes, genomic context)
    +
Machine Learning classifiers
    ↓
Identifies novel variants 
```

**Key Advantages**:
- ✅ Combines known features with ML
- ✅ Identifies novel variants of known virus families
- ✅ Classifies virus types

#### 3. viralFlye - Function-Based Discovery 

```
Pfam protein domain annotation
    ↓
Analyzes domain combination patterns
    ↓
Identifies viruses by FUNCTION (not sequence) 
```

**Key Advantages**:
- ✅ **Function-based** - discovers viruses with novel sequences but known protein functions
- ✅ **Distinguishes** viruses from bacteria/plasmids
- ✅ **Validates** ML predictions with functional evidence

### Three-Tool Synergy for Maximum Discovery Power 

```
Step 1: DeepVirFinder (ML) 
  → Casts wide net, discovers ~200-500 candidates
  → Includes many novel viruses
  
Step 2: VirSorter2 (Hybrid ML) 
  → Validates with features + ML
  → Adds ~30-50 viruses
  
Step 3: viralFlye (Function) 
  → Validates with protein domains
  → Confirms ~28 viruses/candidates
  
Result: Cross-Validation
  → 3-tool consensus: Novel viruses with STRONGEST evidence 
  → 2-tool consensus: Novel viruses with HIGH confidence 
  → ML finds → Function validates → Publication-ready! 
```

**Real-World Impact**:
- Discover viruses in **under-explored environments**
- Identify viruses from **unknown host organisms**
- Find viruses with **no cultured representatives**
- Characterize **dark matter viruses** in metagenomes

---

## Installation Requirements

### System Requirements

- **OS**: Linux/Unix (Recommended: CentOS 7+, Ubuntu 18.04+)
- **Scheduler**: SLURM (optional, also supports local execution)
- **Memory**: 
  - Short-read mode: ≥ 64 GB (SPAdes may require 512 GB)
  - Long-read mode: ≥ 128 GB
  - Long-read + viralFlye: ≥ 256 GB (recommended)
- **Storage**: Reserve 50-200 GB depending on data size

### Core Software

| Software | Minimum Version | Installation |
|----------|----------------|--------------|
| Nextflow | 22.10.0 | conda/official |
| Apptainer/Singularity | 1.0.0 | System admin |
| Conda/Mamba | Latest | Miniconda |
| Python | 3.7+ | conda |
| Java | 11+ | conda/system |

### Detailed Installation Steps

#### 1. Install Nextflow

```bash
# Method 1: conda (recommended)
conda install -c bioconda nextflow

# Method 2: Direct download
curl -fsSL https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

#### 2. Create Conda Environments

```bash
# Environment 1: Main workflow (VirSorter2)
conda create -n nextflow_env -c bioconda -c conda-forge \
    nextflow virsorter2 python=3.9 pandas numpy

# Environment 2: DeepVirFinder
conda create -n dvf -c bioconda -c conda-forge \
    python=3.7 scikit-learn=0.22.1 theano keras=2.3.1 h5py
# Then install following DeepVirFinder official docs

# Environment 3: viralFlye (required for long-read mode)
conda create -n viralFlye_env -c bioconda \
    viralflye flye python=3.10 hmmer
```

#### 3. Download Required Databases

```bash
# VirSorter2 database (~12 GB)
conda activate nextflow_env
virsorter setup -d /path/to/virsorter2/db -j 4

# Pfam database (required for viralFlye, ~2 GB)
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm  # Generate index files

# DeepVirFinder models (included with software)
# Installed with DeepVirFinder
```

#### 4. Install DeepVirFinder

```bash
conda activate dvf
cd /path/to/install/
git clone https://github.com/jessieren/DeepVirFinder.git
cd DeepVirFinder
# Configure following official README
```

#### 5. Configure Workflow

Edit `metagenome_assembly_classification.config`:

```groovy
params {
    // Database paths
    virsorter2_db = '/your/path/to/virsorter2/db'
    deepvirfinder_dir = '/your/path/to/DeepVirFinder'
    pfam_db = '/your/path/to/Pfam-A.hmm'
    
    // Environment paths
    viralflye_env = '/your/path/to/.conda/envs/viralFlye_env'
}

// SLURM configuration
process {
    executor = 'slurm'
    queue = 'your_partition_name'  // Modify to your queue name
}

// Apptainer cache
apptainer {
    cacheDir = '/your/scratch/singularity_cache'
}
```

---

## Input Data Format

### Short-read Data (Illumina)

**File**: `samplesheet.csv` or `samplesheet_short.csv`

```csv
sample,fastq_1,fastq_2
sample1,/absolute/path/to/sample1_R1.fastq.gz,/absolute/path/to/sample1_R2.fastq.gz
sample2,/absolute/path/to/sample2_R1.fastq.gz,/absolute/path/to/sample2_R2.fastq.gz
```

**Requirements**:
- Column names: `sample`, `fastq_1`, `fastq_2`
- Paths: Absolute paths or relative to working directory
- Format: FASTQ (can be gzip compressed)

### Long-read Data (Nanopore/PacBio)

**File**: `samplesheet_long.csv`

```csv
sample,fastq_long
sample1,/absolute/path/to/sample1_nanopore.fastq.gz
sample2,/absolute/path/to/sample2_pacbio.fastq.gz
```

**Requirements**:
- Column names: `sample`, `fastq_long`
- Single-end data (one file)
- Specify platform: `--longread_platform nano` or `pacbio`

---

## Usage

### Method 1: Using SLURM Script (Recommended) 

#### Short-read Data

```bash
# 1. Edit script (if needed)
vim run_metagenome_assembly_classification_shortread.sh

# 2. Submit job
sbatch run_metagenome_assembly_classification_shortread.sh
```

#### Long-read Data + Three-Tool Analysis

```bash
# 1. Edit script, ensure viralFlye is enabled
vim run_metagenome_assembly_classification_longread.sh

# Ensure this line is set to:
ENABLE_VIRALFLYE="true"

# 2. Submit job
sbatch run_metagenome_assembly_classification_longread.sh
```

### Method 2: Nextflow Command Line

#### Basic Long-read Analysis

```bash
nextflow run metagenome_assembly_classification_workflow.nf \
    -c metagenome_assembly_classification.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    --longread true \
    --longread_platform nano
```

#### Complete Three-Tool Analysis (Recommended) 

```bash
nextflow run metagenome_assembly_classification_workflow.nf \
    -c metagenome_assembly_classification.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    --longread true \
    --longread_platform nano \
    --enable_viralflye true \
    --pfam_db /path/to/Pfam-A.hmm \
    --viralflye_env /path/to/.conda/envs/viralFlye_env \
    --viralflye_min_length 500 \
    --deepvirfinder_pvalue 0.05
```

### Method 3: Resume from Failure

```bash
# Use -resume to continue (if only task failed, code unchanged)
nextflow run metagenome_assembly_classification_workflow.nf \
    -c metagenome_assembly_classification.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    -resume
```

⚠️ **Important**: If workflow code was modified, **must clean cache**:
```bash
rm -rf work/ .nextflow*
# Then re-run without -resume
```

---

## Parameter Configuration

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Input sample sheet (CSV) | `samplesheet_long.csv` |
| `--outdir` | Output directory | `results_long` |
| `--virsorter2_db` | VirSorter2 database path | `/data/virsorter2/db` |

### Long-read Mode Parameters

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `--longread` | Enable long-read mode | `false` | `true` |
| `--longread_platform` | Platform type | `nano` | `nano`/`pacbio` |
| `--skip_longread_qc` | Skip long-read QC | `true` | `true` |

### Three-Tool Analysis Parameters (Long-read Only) 

**Note**: viralFlye **only supports long-read mode**

| Parameter | Description | Default | Optimized  |
|-----------|-------------|---------|-------------|
| `--enable_viralflye` | Enable viralFlye (long-read) | `false` | `true` |
| `--viralflye_min_length` | Min viral length (bp) | `1000` | **`500`**  |
| `--viralflye_completeness` | Completeness cutoff | `0.5` | **`0.3`**  |
| `--pfam_db` | Pfam database path | - | Required |
| `--viralflye_env` | viralFlye environment path | - | Required |

### Viral Identification Threshold Parameters (Short+Long)

| Parameter | Description | Default | Optimized  | Mode Support |
|-----------|-------------|---------|-------------|--------------|
| `--virsorter2_min_score` | VirSorter2 min score | `0.5` | `0.5` | Short+Long |
| `--virsorter2_min_length` | VirSorter2 min length | `1000` | `1000` | Short+Long |
| `--deepvirfinder_pvalue` | DeepVirFinder p-value | `0.05` | **`0.05`**  | Short+Long |
| `--deepvirfinder_min_length` | DVF min length | `1000` | `1000` | Short+Long |

**Optimization Notes** (Long-read mode) :
- `viralflye_min_length = 500`: Lower length threshold
- **Automatically includes "Uncertain - viral or bacterial" sequences** 
  - viralFlye reports: 2 → **~28** (2 confirmed + ~26 candidates)
  - High-confidence viruses selected through three-tool validation
- `deepvirfinder_pvalue = 0.05`: Higher sensitivity (more viruses, 200-500)

### QC Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_fastp` | Skip QC | `false` |
| `--save_clean_reads` | Save filtered reads | `true` |
| `--fastp_qualified_quality` | Min quality score | `20` |
| `--fastp_min_length` | Min read length | `50` |

---

## Output Structure

### Long-read Mode + Three-Tool Analysis (Complete Output) 

```
results_long/
│
├── 【Assembly Results】
├── assembly_metaflye/                    # metaFlye contigs (FASTA)
│   └── sample_metaflye_contigs.fa
│
├── metaflye_full_output/                 # metaFlye complete output 
│   └── sample_flye_output/
│       ├── assembly.fasta                # Raw assembly
│       ├── assembly_info.txt             # Assembly stats (depth, circular markers) 
│       ├── assembly_graph.gfa            # Assembly graph (Bandage visualization) 
│       ├── assembly_graph.gv
│       ├── params.json                   # Run parameters
│       ├── flye.log                      # Run log
│       └── [00-assembly/, 10-consensus/, ...]
│
├── 【Viral Identification - Three Parallel Methods】
├── virsorter2_metaflye/                  # VirSorter2 identification
│   ├── sample_vs2_final-viral-score.tsv      # Score table
│   └── sample_vs2_final-viral-combined.fa    # Viral sequences
│
├── deepvirfinder_metaflye/               # DeepVirFinder identification
│   └── sample_dvf_output.txt                 # Prediction results (score, p-value)
│
├── viralflye_results/                    # viralFlye identification (Pfam validated) 
│   ├── sample_viralflye_contigs.fa           # Viral sequences (~28)
│   └── sample_viralflye_summary.csv          # Identification summary
│
├── viralflye_full_output/                # viralFlye complete output
│   └── sample_viralflye_output/
│       ├── vv_circulars/                     # Circular viruses
│       │   ├── circulars_result_table.csv    # Classification results
│       │   └── circulars.fasta               # Viral sequences
│       ├── vv_linears/                       # Linear viruses
│       │   ├── linears_result_table.csv
│       │   └── linears.fasta
│       ├── vc_circulars/                     # viralComplete completeness assessment
│       └── vc_linears/
│
└── 【Three-Tool Comprehensive Comparison】
    └── three_tools_comparison/           # Comprehensive analysis
        ├── sample_three_tools_comparison.txt          # Text report 
        ├── sample_three_tools_comparison.csv          # Detailed data (Excel-ready)
        └── sample_high_confidence_viruses.txt         # High-confidence virus list
```

### Short-read Mode Output

```
results/
├── fastp/                        # QC reports
├── clean_reads/                  # Filtered reads
├── assembly_megahit/             # MEGAHIT assembly
├── assembly_spades/              # metaSPAdes assembly
├── virsorter2_megahit/           # VirSorter2 (MEGAHIT)
├── virsorter2_spades/            # VirSorter2 (SPAdes)
├── deepvirfinder_megahit/        # DeepVirFinder (MEGAHIT)
├── deepvirfinder_spades/         # DeepVirFinder (SPAdes)
├── merged_viral_reports_megahit/ # Merged reports (MEGAHIT)
├── merged_viral_reports_spades/  # Merged reports (SPAdes)
├── assembler_comparison/         # Assembler comparison
└── abundance/                    # Viral abundance analysis (RPM & RPKM) ⭐⭐⭐
    ├── megahit/                 # MEGAHIT viral abundance
    │   ├── *_megahit_abundance.csv
    │   └── *_megahit_abundance_summary.txt
    └── spades/                  # SPAdes viral abundance
        ├── *_spades_abundance.csv
        └── *_spades_abundance_summary.txt
```

---

## Results Interpretation

### Three-Tool Comprehensive Comparison Report (Most Important) 

**File**: `results_long/three_tools_comparison/sample_three_tools_comparison.txt`

**Example Content**:

```
================================================================================
Three-Tool Viral Identification Comprehensive Comparison Report
VirSorter2 + DeepVirFinder + viralFlye (Scheme A: Parallel Independent Analysis)
Sample: llnl_66d1047e
================================================================================

[Overall Statistics]
------------------------------------------------------------------------------------
VirSorter2 identified:        48 viruses
DeepVirFinder identified:     187 viruses (p<0.05)
viralFlye identified:         28 viruses (Pfam validated)
  - Virus (confirmed): 2
  - Uncertain (candidates): 26
Total viruses (deduplicated):  220

[Tool Intersection Analysis]
------------------------------------------------------------------------------------
VirSorter2 ∩ DeepVirFinder:     25
VirSorter2 ∩ viralFlye:         12
DeepVirFinder ∩ viralFlye:      15
3-tool consensus :         10 (Highest confidence)

[Confidence Stratification]
------------------------------------------------------------------------------------
Highest confidence (3-tool):   10 ⭐⭐⭐
Medium confidence (2-tool):    45 ⭐⭐
Low confidence (1-tool):       165 ⭐

[Recommended Analysis Strategy]
------------------------------------------------------------------------------------
1. Prioritize 3-tool consensus viruses (most reliable)
2. 2-tool consensus viruses as secondary priority
3. viralFlye-only viruses (Pfam validated, high specificity)
4. Single-tool viruses for exploratory analysis
```

### High-Confidence Virus List

**File**: `results_long/three_tools_comparison/sample_high_confidence_viruses.txt`

**Format**:

```
# High-confidence viral sequences (prioritized)
# Sample: llnl_66d1047e

# 3-tool consensus (highest confidence): 10
contig_1085    3-tool-consensus
contig_1192    3-tool-consensus
contig_2345    3-tool-consensus
...

# 2-tool consensus: 45
contig_234     2-tool-consensus    VirSorter2+DeepVirFinder
contig_567     2-tool-consensus    VirSorter2+viralFlye
contig_890     2-tool-consensus    DeepVirFinder+viralFlye
...
```

### Detailed Data Table (CSV)

**File**: `results_long/three_tools_comparison/sample_three_tools_comparison.csv`

**Column Descriptions**:

| Column | Description | Example |
|--------|-------------|---------|
| sequence_name | Contig name | contig_1085 |
| identified_by | Tool combination | VirSorter2+DeepVirFinder+viralFlye |
| consensus_count | Consensus level (1-3) | 3 |
| vs2_score | VirSorter2 score | 0.95 |
| vs2_group | VirSorter2 virus type | dsDNAphage |
| dvf_score | DeepVirFinder score | 0.87 |
| dvf_pvalue | DeepVirFinder p-value | 0.0023 |
| vf_score | viralFlye score | 33.68 |
| vf_type | viralFlye type | Circular |
| vf_prediction | viralFlye Pfam prediction | Virus/Uncertain |

**Sorted by `consensus_count` descending - top entries are most reliable!**

---

### Viral Abundance Analysis (RPM & RPKM) ⭐⭐⭐

**Location**: `results/abundance/`

The workflow calculates viral abundance metrics for all identified viruses using two normalized metrics:

#### 1. RPM (Reads Per Million)
- **Formula**: `(reads mapping to contig / total reads) × 1,000,000`
- **Use case**: Compare relative abundance across samples
- **Normalized by**: Total read count only

#### 2. RPKM (Reads Per Kilobase per Million)
- **Formula**: `(reads mapping to contig / contig length in kb) / (total reads / 1,000,000)`
- **Use case**: Compare abundance of contigs with different lengths
- **Normalized by**: Both contig length and total read count
- **Best for**: Identifying the most abundant viral species

#### Output Files

**Short-read Mode**:
```
results_short/abundance/
├── megahit/
│   ├── sample_megahit_abundance.csv           # Detailed abundance data
│   └── sample_megahit_abundance_summary.txt   # Top 10 most abundant viruses
└── spades/
    ├── sample_spades_abundance.csv
    └── sample_spades_abundance_summary.txt
```

**Long-read Mode**:
```
results_long/abundance/
├── metaflye/
│   ├── sample_metaflye_abundance.csv
│   └── sample_metaflye_abundance_summary.txt
└── viralflye/
    ├── sample_viralflye_abundance.csv
    └── sample_viralflye_abundance_summary.txt
```

#### CSV File Format

| Column | Description | Example |
|--------|-------------|---------|
| contig_name | Contig identifier | contig_1085 |
| length | Contig length (bp) | 15432 |
| mapped_reads | Number of reads mapped | 1250 |
| rpm | Reads Per Million | 125.5 |
| rpkm | Reads Per Kilobase per Million | 813.2 |

#### Summary Report

The summary report (`*_abundance_summary.txt`) contains:
- Top 10 most abundant viruses (sorted by RPKM)
- Total viral abundance statistics
- Recommended analysis steps

**Example**:
```bash
cat results_short/abundance/megahit/sample_megahit_abundance_summary.txt
```

#### Interpretation Tips

1. **Use RPKM for comparison**: RPKM accounts for contig length, making it better for comparing different viruses
2. **Compare across assemblers**: Check if the same viruses are abundant in both MEGAHIT and SPAdes results (cross-validation)
3. **Focus on high-confidence viruses**: Combine abundance data with consensus analysis (prioritize 3-tool consensus viruses with high RPKM)
4. **Cross-sample comparison**: Compare RPM values across different samples (if you have multiple samples)

For detailed documentation, see [ABUNDANCE_CALCULATION_README.md](./ABUNDANCE_CALCULATION_README.md)

---

## Analysis Strategy

### Recommended Analysis Workflow

#### 1. View Three-Tool Comparison Summary

```bash
cat results_long/three_tools_comparison/sample_three_tools_comparison.txt
```

#### 2. Extract 3-Tool Consensus Viruses (Most Reliable) 

```bash
# Extract contig IDs
awk -F'\t' '$2=="3-tool-consensus" {print $1}' \
    results_long/three_tools_comparison/sample_high_confidence_viruses.txt \
    > consensus_3_ids.txt

# Extract sequences
seqkit grep -f consensus_3_ids.txt \
    results_long/assembly_metaflye/sample_metaflye_contigs.fa \
    > consensus_3_viruses.fa

# These are the most reliable viruses, prioritize for:
# - Genome annotation (Prokka, DRAM-v)
# - Functional analysis
# - Phylogenetic analysis
# - Host prediction
```

#### 3. Analyze viralFlye-Identified Viruses (Pfam Validated) 

```bash
# Directly use viralFlye output (~28 sequences)
cp results_long/viralflye_results/sample_viralflye_contigs.fa \
   viralflye_verified_viruses.fa

# These viruses:
# - Pfam protein domain validated
# - High specificity
# - May be complete viral genomes
```

#### 4. Analyze 2-Tool Consensus (Medium Confidence) 

```bash
# Extract 2-tool consensus
awk -F'\t' '$2=="2-tool-consensus" {print $1}' \
    results_long/three_tools_comparison/sample_high_confidence_viruses.txt \
    > consensus_2_ids.txt

seqkit grep -f consensus_2_ids.txt \
    results_long/assembly_metaflye/sample_metaflye_contigs.fa \
    > consensus_2_viruses.fa
```

#### 5. Analyze Viral Abundance ⭐⭐⭐

```bash
# View abundance summary (top 10 most abundant viruses)
cat results_long/abundance/metaflye/sample_metaflye_abundance_summary.txt

# Combine with consensus analysis (most abundant + high confidence)
# Filter abundance CSV for 3-tool consensus viruses
awk -F',' 'NR==1 || $1 in consensus_ids' \
    consensus_3_ids.txt \
    results_long/abundance/metaflye/sample_metaflye_abundance.csv \
    > consensus_3_abundance.csv

# These viruses are both:
# - High confidence (3-tool consensus)
# - Highly abundant (high RPKM)
# → Best candidates for detailed analysis
```

---

## Key Concepts

### What is viralFlye?

**Correct Understanding** ✅:
- Viral **identification and validation** tool
- Identifies viruses from metaFlye assembly results
- Uses **Pfam database** to validate viral proteins

**Common Misconceptions** ❌:
- Not an assembly tool
- Not a refinement tool
- Does not re-assemble sequences

### Pfam Database Explained

**Pfam ≠ Virus Protein Database**

**Pfam is**:
- Universal protein family database (~20,000 families)
- Contains proteins from **all life forms** (viruses + bacteria + eukaryotes, etc.)
- HMM profile format

**How viralFlye uses Pfam**:
1. Scans all proteins in contigs
2. Identifies Pfam domains
3. Judges based on domain **combination patterns**:
   - Viral feature domains (Phage_portal, TerL, etc.) → Virus ✅
   - Bacterial markers (Ribosomal_*, DNA_gyrase, etc.) → Chromosome
   - Plasmid features (VirB, Tra system, etc.) → Plasmid

**Advantage**: Can identify **novel, unknown viruses** (as long as they have viral characteristic proteins)

### Why Does viralFlye Identify So Few?

**Reasons**:
1. **Strict standards**: Requires clear viral protein features
2. **Pfam classification bottleneck**: Most candidates classified as Chromosome/Plasmid
3. **High specificity**: Low false positive rate

**Results**:
- Small number identified
- But extremely high quality 
- Low false positive rate

**Strategy for Novel Virus Discovery** :
- **DeepVirFinder**: Many candidates (ML pattern recognition) - **Best for novel viruses** 
- **VirSorter2**: Balanced (ML + features) - Identifies novel variants 
- **viralFlye**: Few but refined (function-based) - Validates by protein domains 
- **Combined three**: **Comprehensive novel virus discovery** 
  - ML finds candidates → Pfam validates function → High confidence

---

## Performance and Resources

### Computational Resource Requirements

| Process | CPU | Memory | Time | Notes |
|---------|-----|--------|------|-------|
| metaFlye | 32 | 128 GB | 24-48h | Depends on data size |
| VirSorter2 | 16 | 64 GB | 12-24h | Depends on contig count |
| DeepVirFinder | 8 | 32 GB | 8-12h | GPU can accelerate |
| viralFlye | 32 | 128 GB | 24-36h | Includes Pfam scan |
| Three-tool comparison | 2 | 8 GB | < 1h | Fast |

### Typical Runtime

| Data Size | Short-read | Long-read | Long+Three-tool |
|-----------|-----------|-----------|-----------------|
| 1 GB | 12-24h | 24-36h | 48-60h |
| 5 GB | 24-48h | 48-72h | 72-96h |
| 10 GB | 48-72h | 72-96h | 96-120h |

### Storage Requirements

| Data Size | Intermediate | Final Results | Total |
|-----------|-------------|---------------|-------|
| 1 GB | ~10 GB | ~5 GB | ~15 GB |
| 5 GB | ~50 GB | ~20 GB | ~70 GB |
| 10 GB | ~100 GB | ~40 GB | ~140 GB |

**Tips**:
- metaFlye complete output takes significant space (~40-50x contigs size)
- Can selectively delete intermediate files after analysis

---

## Frequently Asked Questions

### Q1: What's the difference between short-read and long-read modes?

**A**: Two modes have different design strategies

**Short-read mode**:
- Assemblers: MEGAHIT + metaSPAdes (dual assemblers)
- Viral identification: VirSorter2 + DeepVirFinder
- Comparison: **Assembler comparison** (MEGAHIT vs SPAdes)
- Specialty: Cross-validation through dual assemblers
- viralFlye: ❌ **Not supported**

**Long-read mode** :
- Assembler: metaFlye
- Viral identification: VirSorter2 + DeepVirFinder + **viralFlye** (three tools)
- Comparison: **Three-tool comparison** (VS2 vs DVF vs viralFlye)
- Specialty: viralFlye provides **Pfam protein validation**
- viralFlye: ✅ **Supported**

**How to choose**:
- Have Illumina data → Use short-read mode
- Have Nanopore/PacBio data → Use long-read mode (recommend enabling viralFlye)

### Q2: viralFlye Only Identifies 2 Viruses? 

**A**: Optimized! Now reports **~28 viral candidates** (includes Uncertain sequences)

**How viralFlye works**:
```
metaFlye assembly (1212 contigs)
    ↓
Identify circular/linear sequences: 218
    ↓
Pfam classification (viralVerify)  Key step!
    ↓
├─ Virus (confirmed): 2 
├─ Uncertain - viral or bacterial: ~26 Viral candidates!
├─ Chromosome (bacterial): ~150 (excluded)
└─ Plasmid: ~40 (excluded)
    ↓
viralFlye reports: 2 + 26 = ~28
```

**Key Understanding** :
- **Confirmed viruses** (2): Pfam has clear viral feature domain combinations
- **Viral candidates** (~26): Pfam uncertain, need VS2/DVF validation
- **Three-tool cross-validation** filters high-confidence viruses

**Optimized Results**:
- viralFlye reports: **~28** (2 confirmed + ~26 candidates) 
- 3-tool consensus: **5-15 high-confidence viruses** 
- 2-tool consensus: **15-25** 

### Q3: Error "VIRSORTER2_VIRALFLYE terminated with error"

**A**: Nextflow is using old code cache!

**Solution**:
```bash
# Completely clean cache
rm -rf work/ .nextflow*

# Re-run
sbatch run_metagenome_assembly_classification_longread.sh
```

⚠️ **Important**: After modifying code, **must** clean cache!

### Q4: Does viralFlye Need Flye Installed?

**A**: Yes! viralFlye depends on Flye tool.

```bash
conda activate viralFlye_env
conda install -c bioconda flye
# Or run
bash install_flye_to_viralflye_env.sh
```

### Q5: Do metaFlye and Flye Need Separate Installation?

**A**: No!

metaFlye is Flye's `--meta` mode:
```bash
conda install -c bioconda flye  # Single installation
flye --meta ...  # Use metaFlye mode
```

### Q6: Why Need Complete metaFlye Output? (Long-read Only)

**A**: Multiple purposes

1. **viralFlye required**: Needs complete directory (not just contigs file)
2. **Quality assessment**: assembly_info.txt contains depth, circular markers, etc.
3. **Visualization**: assembly_graph.gfa can be viewed with Bandage
4. **Further analysis**: Can re-run viralFlye or other tools

### Q7: Why Only 2 "Virus" but Report ~28? 

**A**: Includes Pfam-uncertain viral candidate sequences

**Pfam Classification Results**:
- **Virus** (2): Clear viral feature domain combinations 
  - Examples: Phage_portal, TerL_ATPase, Phage_capsid
  - Extremely high confidence, no additional validation needed
  
- **Uncertain - viral or bacterial** (~26): 
  - Pfam cannot determine (too few domains or mixed features)
  - Could be novel viruses, viral fragments, or non-viruses
  - **Require VS2/DVF cross-validation**

**Three-Tool Validation Strategy**:
```
Uncertain sequences + VirSorter2 + DeepVirFinder
  ↓
Identified by 2 or 3 tools → High-confidence virus 
Only viralFlye identifies → Need further validation 
```

**Effect**:
- viralFlye total reports: **~28** (2 confirmed + 26 candidates)
- High-confidence viruses after 3-tool validation: **15-25** 

### Q8: How to Increase Virus Count? (Already Optimized) 

**A**: Workflow automatically optimized

**Optimization Strategy**:
```groovy
// 1. viralFlye: Include Uncertain sequences
viralflye_min_length = 500        // Lower length threshold
// Automatically includes "Uncertain - viral or bacterial"

// 2. DeepVirFinder: Moderate threshold
deepvirfinder_pvalue = 0.05       // Higher sensitivity (more viruses)

// 3. VirSorter2: Standard threshold  
virsorter2_min_score = 0.5
```

**Actual Effects**:
- viralFlye: 2 → **~28**  (includes Uncertain)
- DeepVirFinder: 100 → 200-500 (p=0.05, higher sensitivity)
- **3-tool consensus**: 10-20 high-confidence viruses 

### Q9: What Does Pfam Database Contain?

**A**: Protein families from all life forms (not just viruses)

Pfam contains:
- Viral protein families (~1000+)
- Bacterial protein families (~10000+)
- Eukaryotic protein families (~8000+)
- Total ~20,000 protein families

viralFlye distinguishes viruses, bacteria, and plasmids by analyzing **domain combination patterns**.

### Q10: What's the Relationship Between Three Tools? (Long-read Only)

**A**: Parallel and independent, not sequential!

```
✅ Correct understanding (parallel):
VS2 ∥ DVF ∥ viralFlye → Comprehensive comparison (long-read mode)

❌ Wrong understanding (sequential):
VS2 → DVF → viralFlye → Re-validation
```

Three tools are **independent viral identification methods** (long-read mode only), each with unique characteristics, complementing each other.

### Q11: Can Short-read Mode Use viralFlye?

**A**: ❌ No! viralFlye **only supports long-read data**

**Reasons**:
- viralFlye requires complete metaFlye (Flye) output directory
- metaFlye/Flye only supports long-read data (Nanopore/PacBio)
- MEGAHIT and SPAdes output structures are incompatible

**Short-read Mode Alternatives**:
- ✅ Use dual-tool validation (VirSorter2 + DeepVirFinder)
- ✅ Use dual-assembler consensus (viruses identified by both MEGAHIT and SPAdes)
- ✅ This already provides good cross-validation

**To use viralFlye**:
- Use long-read sequencing data (Nanopore or PacBio)

---

## Advanced Usage

### Extract Viruses by Confidence Level

```bash
# 3-tool consensus (highest confidence)
awk -F'\t' '$2=="3-tool-consensus" {print $1}' \
    results_long/three_tools_comparison/sample_high_confidence_viruses.txt

# viralFlye-only (Pfam validated, high specificity)
awk -F',' '$6=="viralFlye" {print $1}' \
    results_long/three_tools_comparison/sample_three_tools_comparison.csv

# DeepVirFinder high confidence (p < 0.001)
awk 'NR>1 && $4<0.001 {print $1}' \
    results_long/deepvirfinder_metaflye/sample_dvf_output.txt
```

### Visualize metaFlye Assembly Graph

```bash
# Install Bandage
conda install -c bioconda bandage

# Generate image
bandage image \
    results_long/metaflye_full_output/sample_flye_output/assembly_graph.gfa \
    assembly_graph.png \
    --height 3000 --width 3000

# Interactive view
bandage load \
    results_long/metaflye_full_output/sample_flye_output/assembly_graph.gfa
```

### Find Circular Viruses (metaFlye Can Identify)

```bash
# metaFlye marks circular contigs
awk 'NR>1 && $4=="Y" {print $1, $2"bp", "coverage="$3}' \
    results_long/metaflye_full_output/sample_flye_output/assembly_info.txt

# Many viruses (especially phages) have circular genomes
```

### Manually Run viralFlye (If Not Enabled in Workflow)

```bash
conda activate viralFlye_env

viralFlye.py \
    --dir results_long/metaflye_full_output/sample_flye_output \
    --hmm /path/to/Pfam-A.hmm \
    --reads /path/to/original_reads.fastq.gz \
    --outdir manual_viralflye_run \
    --threads 32 \
    --min_viral_length 500 \
    --completeness 0.3
```

### Statistics Tool Performance

```bash
python3 << 'EOF'
import pandas as pd

df = pd.read_csv('results_long/three_tools_comparison/sample_three_tools_comparison.csv')

print("="*60)
print("Three-Tool Identification Statistics")
print("="*60)
print(f"VirSorter2:      {df['vs2_score'].notna().sum():>4} viruses")
print(f"DeepVirFinder:   {df['dvf_score'].notna().sum():>4} viruses")
print(f"viralFlye:       {df['vf_score'].notna().sum():>4} viruses")
print()

print("Consensus Analysis:")
for count in [3, 2, 1]:
    n = len(df[df['consensus_count'] == count])
    stars = '⭐' * count
    print(f"{count}-tool consensus:    {n:>4} viruses {stars}")
print("="*60)
EOF
```

---

## Best Practices

### Parameter Recommendations (Different Research Goals)

#### High-Quality Viral Genome Research

```groovy
// Goal: Few high-quality complete viral genomes
viralflye_min_length = 3000       // Longer
viralflye_completeness = 0.5      // Strict
deepvirfinder_pvalue = 0.001      // Very strict
virsorter2_min_score = 0.7        // High score

// Expected: 5-15 extremely reliable viruses
```

#### Balanced Research (Recommended) 

```groovy
// Goal: Balance quantity and quality
viralflye_min_length = 500        //  Current config
viralflye_completeness = 0.3      //  Current config
deepvirfinder_pvalue = 0.05       //  Current config
virsorter2_min_score = 0.5

// Expected: 50-100 reliable viruses
```

#### Viral Diversity Survey

```groovy
// Goal: Discover as many viruses as possible
viralflye_min_length = 300        // Short
viralflye_completeness = 0.2      // Loose
deepvirfinder_pvalue = 0.1        // Very loose
virsorter2_min_score = 0.3

// Expected: 100-300 viruses (including potential)
```

### Workflow Process Recommendations

1. **First run**: Use default/recommended parameters
2. **View results**: Check identification quantity and quality
3. **Adjust parameters**: Optimize based on needs
4. **Re-run**: Clean cache then run
5. **Result validation**: Use three-tool comparison report

---

## Expected Results (ML-Enhanced Novel Virus Discovery) 

### Typical Sample (Long-read Mode + Three Tools)

```
Input: 10 GB Nanopore metagenomic data

Output:
- metaFlye contigs: ~1000
- VirSorter2 identified: ~30-50 viruses (ML + features)
- DeepVirFinder identified: ~200-500 viruses (Deep Learning, p<0.05) 
- viralFlye identified: ~28 viruses/candidates (Pfam function validation) 
  - Virus (Pfam-confirmed): 2 
  - Uncertain (candidates): ~26 

Three-Tool Comprehensive Validation :
- 3-tool consensus: ~10-20  (ML + Function validated)
- 2-tool consensus: ~30-60 
  - ML-discovered + Function-validated candidates
  - Includes novel viruses with high confidence
- Total usable viruses: ~50-100 (stratified by confidence)

Novel Virus Discovery Potential :
- DeepVirFinder ML: Identifies ~50-100 potential novel viruses
- Cross-validated by other tools: ~20-40 high-confidence novel viruses 
- Pfam-validated novel viruses: ~5-15 (function confirmed) 

Runtime: 48-72 hours
```

**Key Improvements** :
- **ML-powered discovery**: DeepVirFinder finds novel viruses without homology requirement 
- **Multi-tool validation**: ML finds → Function validates → High confidence 
- **viralFlye optimization**: Now reports **~28** (includes Uncertain sequences) 
- **Novel virus focus**: ~20-40 high-confidence novel viruses through cross-validation 
- **Confidence stratification**: Clear tiers (Virus vs Uncertain, 1-3 tool consensus)

---


### Utility Scripts

- **`fix_viralflye_results.sh`** - Quick fix viralFlye results
- **`diagnose_viralflye.sh`** - Environment diagnostics
- **`test_viralflye.sh`** - Tool testing
- **`install_flye_to_viralflye_env.sh`** - Flye installation

---

## Workflow Execution Checklist

### Pre-run Checks

- [ ] Conda environments created (nextflow_env, dvf, viralFlye_env)
- [ ] VirSorter2 database downloaded
- [ ] Pfam database downloaded and indexed (hmmpress)
- [ ] Flye installed in viralFlye_env
- [ ] Sample sheet format correct
- [ ] Config file paths updated
- [ ] Sufficient computational resources

### Runtime Monitoring

```bash
# Check job status
squeue -u $USER

# Real-time log view
tail -f Viral_Classification_LongRead_*.out

# View Nextflow log
tail -f .nextflow.log
```

### Post-run Validation

```bash
# Check key outputs
ls -lh results_long/three_tools_comparison/
ls -lh results_long/viralflye_results/
ls -lh results_long/metaflye_full_output/

# Check identification counts
echo "VirSorter2:"
tail -n +2 results_long/virsorter2_metaflye/*_viral-score.tsv | wc -l

echo "DeepVirFinder (p<0.05):"
awk 'NR>1 && $3<0.05 {count++} END {print count}' \
    results_long/deepvirfinder_metaflye/*_dvf_output.txt

echo "viralFlye:"
grep -c ">" results_long/viralflye_results/*_viralflye_contigs.fa

# View three-tool comparison
cat results_long/three_tools_comparison/*_comparison.txt
```

---

## Troubleshooting

### Diagnostic Tools

```bash
# Diagnose viralFlye environment
bash diagnose_viralflye.sh

# Test viralFlye.py
bash test_viralflye.sh

# Fix viralFlye result extraction
bash fix_viralflye_results.sh
```

### Common Errors and Solutions

#### 1. "VIRSORTER2_VIRALFLYE terminated"

**Cause**: Using old code cache

**Solution**:
```bash
rm -rf work/ .nextflow*
sbatch run_metagenome_assembly_classification_longread.sh
```

#### 2. "viralFlye.py not found"

**Cause**: viralFlye environment configuration issue

**Solution**:
```bash
conda activate viralFlye_env
conda install -c bioconda viralflye flye
```

#### 3. "Pfam database not found"

**Cause**: Pfam path incorrect or not indexed

**Solution**:
```bash
# Check files
ls -lh /path/to/Pfam-A.hmm*

# Re-index
hmmpress Pfam-A.hmm
```

#### 4. VirSorter2 Database Error

**Solution**:
```bash
conda activate nextflow_env
virsorter setup -d /path/to/virsorter2/db -j 4
```

#### 5. Out of Memory

**Solution**: Edit config file
```groovy
withName: 'METAFLYE_ASSEMBLY' {
    memory = '256 GB'  // Increase memory
}
```

---

## Performance Optimization

### 1. Use Fast Storage

```bash
# Set work directory to SSD
export NXF_WORK=/fast/local/scratch/work
```

### 2. Limit Parallel Tasks

```bash
nextflow run ... -process.maxForks 4
```

### 3. Use Mamba for Acceleration

```groovy
conda {
    useMamba = true
}
```

### 4. Clean Intermediate Files

```bash
# After completion
cd results_long/metaflye_full_output/sample_flye_output/
# Delete intermediate directories, keep key files
rm -rf 00-assembly/ 10-consensus/ 20-repeat/ 30-contigger/ 40-polishing/
```

---

## Citation

If using this workflow for publication, please cite the following tools:

### Assembly Tools

- **MEGAHIT**: Li et al. (2015) *Bioinformatics* 31:1674-1676. [doi:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)
- **SPAdes**: Bankevich et al. (2012) *J Comput Biol* 19:455-477. [doi:10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021)
- **Flye/metaFlye**: Kolmogorov et al. (2019) *Nat Biotechnol* 37:540-546. [doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)

### Viral Identification Tools

- **VirSorter2**: Guo et al. (2021) *Microbiome* 9:37. [doi:10.1186/s40168-020-00990-y](https://doi.org/10.1186/s40168-020-00990-y)
- **DeepVirFinder**: Ren et al. (2020) *Quantitative Biology* 8:64-70. [doi:10.1007/s40484-019-0187-4](https://doi.org/10.1007/s40484-019-0187-4)
- **viralFlye/viralVerify**: Antipov et al. (2020) *Bioinformatics* 36:4584-4586. [doi:10.1093/bioinformatics/btaa490](https://doi.org/10.1093/bioinformatics/btaa490)

### Databases and Tools

- **Pfam**: Mistry et al. (2021) *Nucleic Acids Res* 49:D412-D419. [doi:10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913)
- **fastp**: Chen et al. (2018) *Bioinformatics* 34:i884-i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- **Nextflow**: Di Tommaso et al. (2017) *Nat Biotechnol* 35:316-319. [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

---

## Version Information

### Current Version: 5.2.1

**Release Date**: 2025-11-04

**Major Improvements**:
- ✅ Three-tool parallel independent analysis architecture (Scheme A)
- ✅ Remove redundant validation of viralFlye results
- ✅ Add COMPARE_THREE_VIRAL_TOOLS comprehensive comparison
- ✅ **viralFlye viral identification optimization** 
  - Include "Uncertain - viral or bacterial" sequences (~26)
  - viralFlye reports increased from 2 to **~28**
  - High-confidence viruses screened through three-tool cross-validation
- ✅ Fix viralFlye result extraction logic (correctly extract circulars.fasta and linears.fasta)
- ✅ Always save complete metaFlye output (assembly_info.txt, assembly_graph.gfa, etc.)
- ✅ Optimized parameter configuration:
  - `viralflye_min_length = 500` (reduced from 1000)
  - `viralflye_completeness = 0.3` (reduced from 0.5, new parameter)
  - `deepvirfinder_pvalue = 0.05` (balance sensitivity)
- ✅ Performance optimization: Save ~33% computation time

### Version History

- **5.2.1** (2025-11-14): Viral abundance calculation (RPM & RPKM) + viralFlye optimization (Uncertain sequences) + p-value adjustment
- **5.2.0** (2025-11-03): Three-tool parallel analysis + parameter optimization
- **5.1.0** (2025-11-02): Add viralFlye support
- **5.0.0** (2025-11-01): Long-read mode support
- **4.0.0**: Dual-assembler parallel analysis
- **3.0.0**: Initial version

---

## Feature Highlights (v5.2.1)

### 1. Three-Tool Parallel Analysis

**Before**: Cascading validation with redundancy
```
VS2 → DVF → viralFlye → VS2 → DVF (redundant)
```

**Now** : Parallel independent, no redundancy
```
VS2 ∥ DVF ∥ viralFlye → Comprehensive comparison
```

### 2. Confidence Stratification

Automatically generate high-confidence virus list:
- 3-tool consensus 
- 2-tool consensus 
- Single-tool identification 

### 3. viralFlye Optimization 

Significantly increase viral identification:
- **Include "Uncertain - viral or bacterial" sequences**
- viralFlye: From 2 to **~28 viruses/candidates**
- High-confidence viruses (3-tool validated): **15-25**

### 4. Optimized Parameters

Improve virus identification count:
- viralFlye: ~28 (includes Uncertain) 
- DeepVirFinder: 200-500 (p=0.05, higher sensitivity)

### 5. Complete Assembly Information

Save complete metaFlye output:
- assembly_info.txt (depth, circular markers, etc.)
- assembly_graph.gfa (assembly graph, visualizable)
- Support for any subsequent analysis

---

## Important Notes

### Must Clean Cache After Code Modification!

```bash
# After each workflow code modification:
rm -rf work/ .nextflow*

# Then re-run, don't use -resume
```

### viralFlye Environment Must Have Flye Installed

```bash
conda activate viralFlye_env
conda install -c bioconda flye

# Or run installation script
bash install_flye_to_viralflye_env.sh
```

### Pfam is a Universal Database

Pfam contains protein families from **all life forms**, not a virus-specific database. viralFlye identifies viruses by analyzing protein domain combination patterns.

---

### Container Configuration

The workflow uses:
- **Apptainer/Singularity** for MEGAHIT and SPAdes
- **Conda** for fastp, VirSorter2, DeepVirFinder, Flye, minimap2, samtools, seqkit

---

## License

This project is licensed under the MIT License.
Documentation and figures are released under CC BY 4.0.

---

## Acknowledgments

Thanks to the development teams of the following tools:

- **VirSorter2 Team** (DOE Joint Genome Institute)
- **DeepVirFinder Team** (Ren Lab, USC)
- **viralFlye/viralVerify Team** (Antipov Lab, St. Petersburg University)
- **Flye Team** (Kolmogorov Lab, NCI)
- **MEGAHIT and SPAdes Teams**
- **Nextflow Team** (Seqera Labs)
- **Pfam Team** (EMBL-EBI)

---

## Contact

- Email: sihua.peng@uga.edu, Workflow code programmer  
- Email: justin.bahl@uga.edu, Project supervisor  
- Suggestion: [Click here!](https://github.com/pengsihua2023/rvdb-viral-metagenome-nf/issues/new)

---

## Related Resources

### Official Documentation

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [VirSorter2 GitHub](https://github.com/jiarong/VirSorter2)
- [DeepVirFinder GitHub](https://github.com/jessieren/DeepVirFinder)
- [viralFlye GitHub](https://github.com/Dmitry-Antipov/viralFlye)
- [Flye GitHub](https://github.com/fenderglass/Flye)
- [Pfam Database](http://pfam.xfam.org/)

### Related Tools

- **CheckV**: Viral quality assessment
- **VIBRANT**: Viral identification and annotation
- **VirFinder**: Another viral identification tool
- **DIAMOND**: Fast sequence alignment (for functional annotation)
- **Bandage**: Assembly graph visualization
- **seqkit**: Sequence processing tool

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────┐
│                    Input Data                            │
│  Illumina (paired-end) or Nanopore/PacBio (single-end)  │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│              QC (Optional)                               │
│  Short: fastp  │  Long: Filtlong (reserved)            │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│                  Assembly                                │
│  Short: MEGAHIT ∥ SPAdes  │  Long: metaFlye            │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│      Viral Identification (Three Parallel Methods)      │
│  ┌───────────┬────────────────┬──────────────────┐     │
│  │VirSorter2 │ DeepVirFinder  │   viralFlye      │     │
│  │Viral DB   │ Deep Learning  │ Pfam Validation  │     │
│  └─────┬─────┴────────┬───────┴────────┬─────────┘     │
└────────┼──────────────┼────────────────┼───────────────┘
         │              │                │
         └──────────────┼────────────────┘
                        ▼
┌─────────────────────────────────────────────────────────┐
│              Three-Tool Comprehensive Comparison         │
│  - Calculate intersections and unions                    │
│  - Stratify by consensus level (1-3 tools)              │
│  - Generate high-confidence virus list                   │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│                 Final Results                            │
│  ⭐⭐⭐ 3-tool consensus: Highest confidence              │
│  ⭐⭐ 2-tool consensus: Medium confidence                │
│  ⭐ Single-tool: Exploratory analysis                   │
└─────────────────────────────────────────────────────────┘
```

---

## Recommended Reading

### Viral Metagenomics

- Paez-Espino et al. (2016) *Nature* 536:425-430.
- Roux et al. (2019) *Nat Microbiol* 4:1895-1906.

### Methodology

- Nayfach et al. (2021) *Nat Biotechnol* 39:103-111. (CheckV)
- Kieft et al. (2020) *PeerJ* 8:e9439. (VIBRANT)

---

**Last Updated**: 2025-11-04  
**Workflow Version**: 5.2.1  
**Design Philosophy**: Three-tool parallel independent analysis, comprehensive comparison, confidence stratification

**Happy viral hunting!** 

---

## Get Started Now

```bash
# 1. Clean cache (if previously run)
rm -rf work/ .nextflow*

# 2. Run workflow
sbatch run_metagenome_assembly_classification_longread.sh

# 3. Wait for completion, view results
cat results_long/three_tools_comparison/*_comparison.txt
```


