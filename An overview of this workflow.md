# AI-Powered Viral Discovery from Metagenomic Data

---

## Slide 1: Title Slide
**Title:** AI-Powered Viral Discovery from Metagenomic Data  
**Subtitle:** Machine Learning & Deep Learning Pipeline for Novel Virus Discovery  
**Presenter:** Sihua Peng  
**Date:** 10/31/2025  
**Institution:** University of Georgia

---

## Slide 2: Research Background
### The Challenge of Viral Discovery in Metagenomes

- **Hidden Viral Diversity**: Metagenomic samples contain vast uncultured viral diversity
- **Low Abundance**: Many viruses exist at low abundance, making detection challenging
- **Novel Viruses**: Unknown viruses lack reference sequences for traditional methods
- **Complex Communities**: Mixed DNA/RNA viruses in environmental samples
- **Traditional Limitations**: Database-dependent methods miss novel viruses

**→ Need for AI-powered, database-independent approaches**

---

## Slide 3: Research Objectives
### What We Aim to Achieve

✅ **Novel Virus Discovery**: Identify previously unknown viral sequences  
✅ **Low-Abundance Detection**: Detect viruses at low abundance levels  
✅ **Comprehensive Coverage**: Support both DNA and RNA viruses  
✅ **High Confidence**: Dual validation using multiple AI methods  
✅ **Automated Pipeline**: Streamlined workflow for large-scale analysis  

---

## Slide 4: Workflow Overview
### Integrated Multi-Tool Pipeline

```
Raw Metagenomic Data
        ↓
   [QC & Assembly]
        ↓
[Viral Identification]
    ↙        ↘
VirSorter2   DeepVirFinder
  (ML)        (DL)
    ↘        ↙
   [Integration & Analysis]
        ↓
   Final Results
```

**Two AI-powered tools working in parallel for comprehensive viral identification**

---

## Slide 5: Core Technologies - VirSorter2
### Machine Learning-Based Viral Identification

**Method:**
- Supervised machine learning models
- Rule-based detection system
- Virus-specific feature recognition

**Capabilities:**
- ✅ Novel virus detection (database-independent)
- ✅ Multiple virus types: dsDNA, ssDNA, RNA, NCLDV
- ✅ Viral hallmark gene identification
- ✅ Confidence scoring (0-1 scale)

**Advantage:** Identifies viruses based on sequence features, not just similarity

---

## Slide 6: Core Technologies - DeepVirFinder
### Deep Learning-Based Viral Prediction

**Method:**
- Deep convolutional neural network
- Trained on viral sequence patterns
- Statistical significance testing

**Capabilities:**
- ✅ Pattern recognition in sequences
- ✅ P-value-based significance filtering
- ✅ Database-independent approach
- ✅ Novel virus discovery capability

**Advantage:** Learns complex patterns that may indicate viral sequences

---

## Slide 7: Dual Validation Strategy
### Consensus-Based High-Confidence Detection

```
VirSorter2 Results     DeepVirFinder Results
     ↓                        ↓
     └────────┬───────────────┘
              ↓
      Consensus Sequences
    (Identified by Both Tools)
              ↓
    High-Confidence Viral Sequences
```

**Benefits:**
- Reduced false positives
- Increased confidence in novel virus calls
- Complementary strengths of ML and DL

---

## Slide 8: Assembly Strategy - Short Reads
### Dual Assembler Approach

**MEGAHIT**
- Fast, memory-efficient
- Good for complex communities
- Optimized for metagenomes

**metaSPAdes**
- High-quality assemblies
- Better for complex regions
- Higher resource requirements

**Advantage:** Cross-validation between assemblers increases reliability

---

## Slide 9: Assembly Strategy - Long Reads
### Long-Read Assembly for Better Resolution

**metaFlye**
- Designed for long-read data (Nanopore/PacBio)
- Better for low-abundance viruses
- Resolves repetitive regions
- Full-length viral genome recovery

**Platform Support:**
- Nanopore (MinION, PromethION)
- PacBio (Sequel, Sequel II)

**Optional: viralFlye Refinement**
- Targeted reassembly of viral contigs
- Improved assembly quality

---

## Slide 10: Workflow Architecture
### Complete Analysis Pipeline

**Stage 1: Quality Control**
- fastp filtering (short-reads)
- Optional long-read QC

**Stage 2: Assembly**
- Short-read: MEGAHIT + metaSPAdes (parallel)
- Long-read: metaFlye

**Stage 3: Viral Identification**
- VirSorter2 analysis
- DeepVirFinder analysis
- (Optional: viralFlye refinement for long-reads)

**Stage 4: Integration**
- Result merging
- Consensus identification
- Comprehensive reporting

---

## Slide 11: Key Features
### Comprehensive Viral Discovery Pipeline

✅ **Novel Virus Discovery**: AI-powered, database-independent  
✅ **Low-Abundance Detection**: Sensitive assembly methods  
✅ **DNA & RNA Support**: Comprehensive virus type coverage  
✅ **Dual Validation**: ML + DL for high confidence  
✅ **Automated Workflow**: Nextflow-based, reproducible  
✅ **Scalable**: HPC cluster optimized (SLURM)  
✅ **Flexible**: Short-read and long-read support  
✅ **Refinement Option**: viralFlye for improved assemblies  

---

## Slide 12: Application Scenarios
### Where This Workflow Shines

**Environmental Metagenomics**
- Soil, water, ocean samples
- Viral diversity exploration

**Human/Animal Microbiomes**
- Gut virome studies
- Pathogen discovery

**Novel Virus Discovery**
- Emerging virus identification
- Viral surveillance

**Low-Abundance Detection**
- Rare virus discovery
- Pathogen monitoring

**Comparative Analysis**
- Multi-sample viral profiling
- Temporal viral dynamics

---

## Slide 13: Detection Capabilities
### What Can Be Detected?

**Virus Types:**
- dsDNA phages
- ssDNA viruses
- RNA viruses
- NCLDV (Nucleocytoplasmic Large DNA Viruses)
- Lavidaviridae (helper viruses)

**Sequence Characteristics:**
- Novel/unknown viruses ✅
- Low-abundance viruses ✅ (assembly-dependent)
- Partial genomes (fragments)
- Complete viral genomes (long-reads)

**Confidence Levels:**
- High: Consensus (both tools)
- Medium: Single tool, high score
- Low: Single tool, low score

---

## Slide 14: Output Results
### Comprehensive Reporting System

**Generated Files:**
1. **Assembly Contigs**: Assembled viral sequences
2. **Viral Scores**: VirSorter2 confidence scores
3. **DL Predictions**: DeepVirFinder scores and p-values
4. **Merged Reports**: Integrated analysis
5. **Consensus Lists**: High-confidence sequences

**Report Types:**
- Text reports (statistics)
- CSV tables (detailed data)
- Consensus sequences (high-confidence list)

---

## Slide 15: Workflow Comparison
### Advantages Over Traditional Methods

| Feature | Traditional (BLAST) | This Workflow |
|---------|-------------------|---------------|
| Novel Viruses | ❌ Limited | ✅ Strong |
| Database Dependency | ❌ Required | ✅ Optional |
| Low Abundance | ⚠️ Limited | ✅ Better |
| DNA Viruses | ✅ Good | ✅ Excellent |
| RNA Viruses | ⚠️ Variable | ✅ Supported |
| Confidence | ⚠️ Single method | ✅ Dual validation |
| Automation | ⚠️ Manual | ✅ Automated |

---

## Slide 16: Technical Advantages
### Why This Approach Works

**AI-Powered Detection**
- Learns patterns, not just matches
- Adapts to new viral sequences
- Reduces false positives

**Dual Assembler Strategy**
- Cross-validation increases reliability
- Different strengths complement each other

**Dual Identification Tools**
- ML + DL = comprehensive coverage
- Consensus = high confidence

**Long-Read Support**
- Better for low-abundance
- Full-length genomes
- Resolves complexity

---

## Slide 17: Usage Example - Short Reads
### Simple Command-Line Execution

```bash
nextflow run workflow.nf \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/db \
    --deepvirfinder_dir /path/to/DeepVirFinder
```

**Samplesheet Format:**
```csv
sample,fastq_1,fastq_2
sample1,reads_R1.fq.gz,reads_R2.fq.gz
```

**Output:**
- MEGAHIT + SPAdes assemblies
- VirSorter2 + DeepVirFinder results
- Merged reports and consensus sequences

---

## Slide 18: Usage Example - Long Reads
### Nanopore/PacBio Data Processing

```bash
nextflow run workflow.nf \
    --input samplesheet_long.csv \
    --outdir results_long \
    --longread true \
    --longread_platform nano \
    --enable_viralflye true
```

**Features:**
- metaFlye assembly
- Optional viralFlye refinement
- Better low-abundance detection
- Full-length viral genomes

---

## Slide 19: Performance & Scalability
### HPC-Optimized for Large-Scale Analysis

**Resource Efficiency:**
- Parallel processing
- SLURM cluster support
- Configurable resource allocation

**Scalability:**
- Multiple samples simultaneously
- Efficient containerization (Apptainer/Singularity)
- Conda environment management

**Reproducibility:**
- Nextflow ensures reproducible results
- Version-controlled workflows
- Containerized tools

---

## Slide 20: Validation & Quality Control
### Ensuring Reliable Results

**Multi-Level Validation:**
1. Dual assembler cross-check
2. Dual AI tool consensus
3. Statistical significance testing
4. Confidence scoring

**Quality Metrics:**
- VirSorter2 scores (0-1)
- DeepVirFinder p-values
- Consensus agreement
- Assembly quality metrics

**Result Filtering:**
- Configurable thresholds
- Length requirements
- Significance filters

---

## Slide 21: Case Study / Example Results
### Real-World Application

**Sample Scenario:**
- Environmental metagenome sample
- [Your specific example if available]

**Results:**
- X viral sequences identified
- Y consensus high-confidence viruses
- Z novel sequences (not in databases)

**Key Findings:**
- [Highlight interesting discoveries]
- [Novel virus examples]
- [Low-abundance detections]

---

## Slide 22: Limitations & Considerations
### Important Caveats

**Assembly Dependencies:**
- Low-abundance detection requires sufficient coverage
- Very short fragments may be missed
- Assembly quality affects detection

**RNA Virus Detection:**
- Primary reliance on VirSorter2
- DeepVirFinder optimized for DNA

**Computational Resources:**
- SPAdes requires high memory (512GB)
- Long-read assembly is time-intensive
- Large samples need adequate storage

**Best Practices:**
- Ensure sufficient sequencing depth
- Use appropriate quality thresholds
- Validate novel discoveries experimentally

---

## Slide 23: Future Directions
### Potential Enhancements

**Planned Improvements:**
- Enhanced RNA virus detection models
- Integration of additional viral prediction tools
- Improved long-read QC integration
- Functional annotation pipeline

**Research Opportunities:**
- Viral-host interaction prediction
- Taxonomic classification
- Phylogenetic analysis
- Comparative genomics

---

## Slide 24: Summary
### Key Takeaways

✅ **Comprehensive Pipeline**: End-to-end viral discovery workflow  
✅ **AI-Powered**: ML + DL for novel virus detection  
✅ **Dual Validation**: High-confidence consensus results  
✅ **Flexible**: Short-read and long-read support  
✅ **Automated**: Reproducible, scalable analysis  
✅ **Broad Coverage**: DNA and RNA viruses  

**Impact:** Enables systematic viral discovery in metagenomic datasets

---







