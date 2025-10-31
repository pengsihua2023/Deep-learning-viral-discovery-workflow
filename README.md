# Metagenome Viral Classification Workflow

A comprehensive Nextflow workflow for identifying viral sequences in metagenomic data using multiple assembly methods and viral identification tools.

## Overview

This workflow integrates state-of-the-art tools for metagenomic viral sequence identification:

1. **Quality Control**: Optional read filtering using `fastp`
2. **Assembly**: Parallel assembly using `MEGAHIT` and `metaSPAdes` (short-reads) or `metaFlye` (long-reads)
3. **Viral Identification**: 
   - `VirSorter2`: Machine learning-based viral sequence identification
   - `DeepVirFinder`: Deep learning-based viral prediction
4. **Result Integration**: Comprehensive merging and comparison of results from multiple tools
5. **Optional Refinement**: Targeted reassembly of viral contigs using `viralFlye` (long-reads only)

## Features

- ✅ **Dual Assembler Support** (short-reads): Parallel assembly with MEGAHIT and metaSPAdes
- ✅ **Long-read Support**: Assembly using metaFlye for PacBio/Nanopore data
- ✅ **Dual Viral Identification**: VirSorter2 + DeepVirFinder for comprehensive coverage
- ✅ **Consensus Detection**: Identifies high-confidence viral sequences agreed upon by both tools
- ✅ **Optional Refinement**: viralFlye for targeted reassembly of viral contigs from long-reads
- ✅ **SLURM Compatible**: Optimized for HPC cluster environments
- ✅ **Container Support**: Uses Apptainer/Singularity for reproducible execution
- ✅ **Conda Integration**: Flexible tool installation via Conda environments

## Requirements

### Software Dependencies

- **Nextflow** (≥22.10.0)
- **Apptainer/Singularity** (for MEGAHIT, SPAdes, Flye containers)
- **Conda/Mamba** (for fastp, VirSorter2, DeepVirFinder, minimap2, samtools, seqkit)
- **Python** (≥3.7, for result merging scripts)
- **SLURM** (for cluster execution, optional)

### Tool Installations

1. **VirSorter2**: Pre-installed in conda environment (`nextflow_env`)
   - Database required: Download from [VirSorter2 GitHub](https://github.com/jiarong/VirSorter2)

2. **DeepVirFinder**: Pre-installed in conda environment (`dvf`)
   - Installation directory must be specified via `--deepvirfinder_dir`

3. **Assembly Tools**:
   - MEGAHIT: Apptainer container
   - SPAdes: Apptainer container  
   - Flye: Conda environment (`bioconda::flye=2.9`)

4. **Additional Tools**:
   - fastp: Conda environment
   - minimap2, samtools, seqkit: Conda environments (for viralFlye branch)

## Installation

1. **Clone or download this repository**

2. **Install Nextflow**:
   ```bash
   conda install -c bioconda nextflow
   # or
   curl -fsSL https://get.nextflow.io | bash
   ```

3. **Set up Conda environments**:
   ```bash
   # Create environment with Nextflow and VirSorter2
   conda create -n nextflow_env -c bioconda nextflow virsorter2
   conda activate nextflow_env
   
   # Create environment with DeepVirFinder
   conda create -n dvf python=3.7
   conda activate dvf
   # Install DeepVirFinder according to its documentation
   ```

4. **Download VirSorter2 database**:
   ```bash
   wget https://zenodo.org/record/4580235/files/virsorter2_v2.2_db.tar.gz
   tar -xzf virsorter2_v2.2_db.tar.gz
   ```

5. **Configure paths** in `metagenome_assembly_classification_en.config`:
   - Update `deepvirfinder_dir` path
   - Adjust resource allocations if needed
   - Update Apptainer cache directory

## Input Data Format

### Short-read Samplesheet (Illumina)

CSV format with header:
```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

**File**: `samplesheet.csv` or `samplesheet_short.csv`

### Long-read Samplesheet (PacBio/Nanopore)

CSV format with header:
```csv
sample,fastq_long
sample1,/path/to/sample1_nanopore.fastq.gz
sample2,/path/to/sample2_pacbio.fastq.gz
```

**File**: `samplesheet_long.csv`

## Usage

### Quick Start

#### Short-read Data (Default)

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder
```

#### Long-read Data (Nanopore)

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    --longread true \
    --longread_platform nano
```

#### Long-read Data (PacBio)

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    --longread true \
    --longread_platform pacbio
```

### Using SLURM Runner Scripts

For convenience, use the provided runner scripts:

#### Short-read Workflow
```bash
sbatch run_metagenome_assembly_classification_shortread_en.sh
```

#### Long-read Workflow
```bash
sbatch run_metagenome_assembly_classification_longread_en.sh
```

**Note**: Edit the runner scripts to update paths and parameters before submission.

### Advanced Usage

#### Enable viralFlye Refinement (Long-reads only)

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --longread true \
    --longread_platform nano \
    --enable_viralflye true \
    --viralflye_min_score 0.6 \
    --viralflye_min_length 1500
```

#### Skip Quality Control

```bash
--skip_fastp true          # Skip fastp QC (short-reads)
--skip_longread_qc true    # Skip long-read QC (default: true)
```

#### Skip Specific Tools

```bash
--skip_virsorter2 true     # Skip VirSorter2 analysis
--skip_deepvirfinder true  # Skip DeepVirFinder analysis
--skip_merge_reports true  # Skip result merging
```

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Input samplesheet (CSV) | `samplesheet.csv` |
| `--outdir` | Output directory | `results` |
| `--virsorter2_db` | VirSorter2 database path | `/path/to/virsorter2/db` |

### Long-read Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--longread` | Enable long-read mode | `false` |
| `--longread_platform` | Platform: `nano` or `pacbio` | `nano` |
| `--skip_longread_qc` | Skip long-read QC | `true` |
| `--enable_viralflye` | Enable viralFlye refinement | `false` |
| `--viralflye_min_score` | Min VS2 score for targets | `0.5` |
| `--viralflye_min_length` | Min contig length for targets | `1000` |

### Quality Control Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_fastp` | Skip fastp QC | `false` |
| `--save_clean_reads` | Save filtered reads | `true` |
| `--fastp_qualified_quality` | Min quality score | `20` |
| `--fastp_unqualified_percent` | Max low-quality bases % | `40` |
| `--fastp_min_length` | Min read length | `50` |

### Viral Identification Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--virsorter2_min_length` | Min contig length | `1000` |
| `--virsorter2_min_score` | Min viral score | `0.5` |
| `--deepvirfinder_min_length` | Min contig length | `1000` |
| `--deepvirfinder_pvalue` | P-value threshold | `0.05` |

### Assembly Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--megahit_memory` | MEGAHIT memory fraction | `0.8` |
| `--megahit_min_contig_len` | Min contig length | `1000` |
| `--spades_meta` | Use metaSPAdes mode | `true` |

## Output Structure

### Short-read Workflow Output

```
results/
├── fastp/                          # Quality control reports
│   ├── sample1_fastp.html
│   └── sample1_fastp.json
├── clean_reads/                    # Filtered reads (if saved)
│   ├── sample1_clean_R1.fastq.gz
│   └── sample1_clean_R2.fastq.gz
├── assembly_megahit/               # MEGAHIT assembly
│   └── sample1_megahit_contigs.fa
├── assembly_spades/                # metaSPAdes assembly
│   └── sample1_spades_contigs.fa
├── virsorter2_megahit/             # VirSorter2 (MEGAHIT)
│   ├── sample1_megahit_vs2_final-viral-score.tsv
│   └── sample1_megahit_vs2_final-viral-combined.fa
├── virsorter2_spades/              # VirSorter2 (SPAdes)
│   ├── sample1_spades_vs2_final-viral-score.tsv
│   └── sample1_spades_vs2_final-viral-combined.fa
├── deepvirfinder_megahit/          # DeepVirFinder (MEGAHIT)
│   └── sample1_megahit_dvf_output.txt
├── deepvirfinder_spades/           # DeepVirFinder (SPAdes)
│   └── sample1_spades_dvf_output.txt
├── merged_viral_reports_megahit/   # Merged reports (MEGAHIT)
│   ├── sample1_megahit_viral_merged_report.txt
│   ├── sample1_megahit_viral_merged_report.csv
│   └── sample1_megahit_viral_consensus.txt
├── merged_viral_reports_spades/    # Merged reports (SPAdes)
│   ├── sample1_spades_viral_merged_report.txt
│   ├── sample1_spades_viral_merged_report.csv
│   └── sample1_spades_viral_consensus.txt
└── assembler_comparison/           # Assembler comparison
    ├── sample1_assembler_comparison.txt
    ├── sample1_assembler_comparison.csv
    └── sample1_consensus_viral_sequences.txt
```

### Long-read Workflow Output

```
results_long/
├── assembly_metaflye/              # metaFlye assembly
│   └── sample1_metaflye_contigs.fa
├── virsorter2_metaflye/            # VirSorter2 results
│   ├── sample1_metaflye_vs2_final-viral-score.tsv
│   └── sample1_metaflye_vs2_final-viral-combined.fa
├── deepvirfinder_metaflye/         # DeepVirFinder results
│   └── sample1_metaflye_dvf_output.txt
└── merged_viral_reports_metaflye/  # Merged reports
    ├── sample1_metaflye_viral_merged_report.txt
    ├── sample1_metaflye_viral_merged_report.csv
    └── sample1_metaflye_viral_consensus.txt
```

### viralFlye Refinement Output (if enabled)

```
results_long/
├── viralflye_targets/              # Selected viral targets
│   ├── sample1_viral_target_ids.txt
│   └── sample1_viral_targets.fa
├── viralflye_reads/                # Extracted viral reads
│   └── sample1_viral_reads.fastq.gz
├── assembly_viralflye/             # Refined assembly
│   └── sample1_viralflye_contigs.fa
├── virsorter2_viralflye/           # Re-annotation (VirSorter2)
│   └── sample1_viralflye_vs2_final-viral-score.tsv
├── deepvirfinder_viralflye/        # Re-annotation (DeepVirFinder)
│   └── sample1_viralflye_dvf_output.txt
└── merged_viral_reports_viralflye/ # Merged reports (refined)
    ├── sample1_viralflye_viral_merged_report.txt
    ├── sample1_viralflye_viral_merged_report.csv
    └── sample1_viralflye_viral_consensus.txt
```

## Result Files

### Merged Report Files

1. **`*_viral_merged_report.txt`**: Comprehensive text report with statistics
   - Number of viral sequences identified by each tool
   - Consensus sequences (identified by both tools)
   - Tool-specific identifications

2. **`*_viral_merged_report.csv`**: Detailed CSV table
   - Sequence name
   - Identification method (Both/VirSorter2_only/DeepVirFinder_only)
   - VirSorter2 scores and groups
   - DeepVirFinder scores and p-values
   - Consensus flag

3. **`*_viral_consensus.txt`**: List of high-confidence viral sequences
   - Sequences identified by both VirSorter2 and DeepVirFinder

### Assembler Comparison Files (Short-reads only)

1. **`*_assembler_comparison.txt`**: Comparison report between MEGAHIT and SPAdes
2. **`*_assembler_comparison.csv`**: Detailed comparison data
3. **`*_consensus_viral_sequences.txt`**: Final high-confidence sequences (both assemblers)

## Workflow Details

### Short-read Workflow

1. **Quality Control** (optional): `fastp` filtering
2. **Assembly**: Parallel execution of MEGAHIT and metaSPAdes
3. **Viral Identification**: VirSorter2 and DeepVirFinder on both assemblies
4. **Result Merging**: Integration of results per assembler
5. **Assembler Comparison**: Cross-assembler consensus identification

### Long-read Workflow

1. **Quality Control** (optional): Long-read filtering (placeholder for Filtlong/NanoFilt)
2. **Assembly**: metaFlye with `--meta` mode
3. **Viral Identification**: VirSorter2 and DeepVirFinder
4. **Result Merging**: Integration of both tool results
5. **Optional Refinement** (viralFlye):
   - Select viral contigs based on VS2 score and length thresholds
   - Extract reads mapping to viral contigs using minimap2
   - Targeted reassembly using Flye
   - Re-annotation with VirSorter2 and DeepVirFinder

## Configuration

### Resource Allocation

Edit `metagenome_assembly_classification_en.config` to adjust:

- CPU cores per process
- Memory allocation
- Time limits
- SLURM queue settings

### Container Configuration

The workflow uses:
- **Apptainer/Singularity** for MEGAHIT and SPAdes
- **Conda** for fastp, VirSorter2, DeepVirFinder, Flye, minimap2, samtools, seqkit

Container cache directory can be configured in the config file.

## Troubleshooting

### Common Issues

1. **Container Pull Failures**
   - Check internet connection
   - Verify Apptainer/Singularity is properly configured
   - Check `pullTimeout` setting in config

2. **Conda Environment Not Found**
   - Verify environment paths in config file
   - Ensure environments are activated correctly
   - Check `CONDA_BASE` path

3. **Database Path Errors**
   - Verify VirSorter2 database path
   - Check DeepVirFinder directory path
   - Ensure paths are accessible from compute nodes

4. **Memory Errors**
   - Increase memory allocation in config
   - Reduce number of parallel processes
   - Consider using higher memory nodes for SPAdes

5. **Empty Input Channels**
   - Verify samplesheet format
   - Check file paths exist and are accessible
   - Ensure no empty rows in samplesheet

### Debugging

- Check `.nextflow.log` for detailed error messages
- Inspect process work directories: `work/<hash>/.command.err`
- Use `-resume` to continue from failed steps
- Enable verbose logging: `nextflow run ... -with-report report.html`

## Citation

If you use this workflow, please cite the following tools:

- **MEGAHIT**: Li et al. (2015) Bioinformatics. [doi:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)
- **SPAdes**: Bankevich et al. (2012) J Comput Biol. [doi:10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021)
- **Flye**: Kolmogorov et al. (2019) Nat Biotechnol. [doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)
- **VirSorter2**: Guo et al. (2021) Microbiome. [doi:10.1186/s40168-020-00990-y](https://doi.org/10.1186/s40168-020-00990-y)
- **DeepVirFinder**: Ren et al. (2020) Bioinformatics. [doi:10.1093/bioinformatics/btaa010](https://doi.org/10.1093/bioinformatics/btaa010)
- **fastp**: Chen et al. (2018) Bioinformatics. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

## Version

- **Workflow Version**: 5.1.0
- **Nextflow**: DSL2

## License

[Specify your license here]

## Contact

For questions or issues, please [open an issue](https://github.com/your-repo/issues) or contact the maintainers.

## Acknowledgments

This workflow integrates multiple open-source tools for metagenomic viral sequence identification. We thank all tool developers for their contributions to the community.

