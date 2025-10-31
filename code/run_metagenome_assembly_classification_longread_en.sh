#!/bin/bash
#SBATCH --job-name=Viral_Classification_LongRead
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Classification_LongRead_%j.out
#SBATCH --error=Viral_Classification_LongRead_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow (Long-Read Mode)"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

# User's conda environment path
USER_CONDA_ENV="/home/sp96859/.conda/envs/nextflow_env"

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "⚠️  Warning: conda info --base failed, trying alternative method..."
    CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
fi

echo "   Conda base: $CONDA_BASE"
echo "   Target env: $USER_CONDA_ENV"

# Initialize conda
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "❌ Cannot find conda.sh at $CONDA_BASE/etc/profile.d/conda.sh"
    exit 1
fi

# Check if user environment exists
if [ ! -d "$USER_CONDA_ENV" ]; then
    echo "❌ Conda environment not found: $USER_CONDA_ENV"
    echo "   Available environments:"
    conda env list
    exit 1
fi

# Activate environment using absolute path
conda activate "$USER_CONDA_ENV"

# Force update PATH
export PATH="$USER_CONDA_ENV/bin:$PATH"

# Set conda-related environment variables
export CONDA_PREFIX="$USER_CONDA_ENV"
export CONDA_DEFAULT_ENV="nextflow_env"

# Get Python version and set PYTHONPATH
PYTHON_VERSION=$("$USER_CONDA_ENV/bin/python" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "3.9")
export PYTHONPATH="$USER_CONDA_ENV/lib/python${PYTHON_VERSION}/site-packages:${PYTHONPATH:-}"

# Verify environment activation
PYTHON_PATH=$(which python)
echo "   After PATH update:"
echo "   - Python path: $PYTHON_PATH"
echo "   - CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"

if [[ "$PYTHON_PATH" == *"$USER_CONDA_ENV"* ]]; then
    echo "✅ Conda environment activated successfully!"
    echo "   Environment: nextflow_env"
    echo "   Python: $PYTHON_PATH"
else
    echo "❌ Failed to activate user conda environment!"
    echo "   Expected Python in: $USER_CONDA_ENV"
    echo "   Actual Python: $PYTHON_PATH"
    exit 1
fi

# Verify tools
echo "🧪 2. Verifying tools..."
echo "✅ Nextflow: $(which nextflow)"

# Check for Apptainer/Singularity
if command -v apptainer &> /dev/null; then
    echo "✅ Apptainer: $(which apptainer)"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(which singularity)"
else
    echo "❌ Apptainer/Singularity not found (required for containers)"
    exit 1
fi

echo ""
echo "ℹ️  Note: Long-read workflow execution environment"
echo "   - Assembly tool: metaFlye (Apptainer container)"
echo "   - Optional QC: Long-read filtering (if enabled)"
echo "   - Optional refinement: viralFlye (if enabled)"
echo "   - Viral identification tools:"
echo "     * VirSorter2: Pre-installed in nextflow_env ✅"
echo "     * DeepVirFinder: Pre-installed in dvf environment ✅"
echo ""

# Set database paths
VIRSORTER2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/virsorter2/db"
DEEPVIRFINDER_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder"

# Verify databases and tools
echo "🗄️ 3. Verifying databases and tools..."

# VirSorter2 database
if [ -d "$VIRSORTER2_DB" ]; then
    echo "✅ VirSorter2 database: $VIRSORTER2_DB"
    echo "   Database size: $(du -sh $VIRSORTER2_DB | cut -f1)"
else
    echo "❌ VirSorter2 database not found: $VIRSORTER2_DB"
    exit 1
fi

# DeepVirFinder installation
if [ -d "$DEEPVIRFINDER_DIR" ] && [ -f "$DEEPVIRFINDER_DIR/dvf.py" ]; then
    echo "✅ DeepVirFinder: $DEEPVIRFINDER_DIR"
else
    echo "❌ DeepVirFinder not found: $DEEPVIRFINDER_DIR"
    exit 1
fi

echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "samplesheet_long.csv" ]; then
    echo "✅ Samplesheet: samplesheet_long.csv"
    echo "📊 Found $(tail -n +2 samplesheet_long.csv | grep -v '^$' | wc -l) samples"
else
    echo "❌ Samplesheet not found: samplesheet_long.csv"
    echo "   Please create samplesheet_long.csv with format:"
    echo "   sample,fastq_long"
    echo "   sample1,/path/to/long_reads.fastq.gz"
    exit 1
fi

# Long-read platform detection (default: nano, can be set via command line or env)
LONGREAD_PLATFORM="${LONGREAD_PLATFORM:-nano}"
ENABLE_VIRALFLYE="${ENABLE_VIRALFLYE:-false}"
VIRALFLYE_MIN_SCORE="${VIRALFLYE_MIN_SCORE:-0.6}"
VIRALFLYE_MIN_LENGTH="${VIRALFLYE_MIN_LENGTH:-1500}"

echo ""
echo "🔧 Long-read configuration:"
echo "   Platform: $LONGREAD_PLATFORM (Nanopore/PacBio)"
echo "   viralFlye refinement: $ENABLE_VIRALFLYE"
if [ "$ENABLE_VIRALFLYE" = "true" ]; then
    echo "   viralFlye min score: $VIRALFLYE_MIN_SCORE"
    echo "   viralFlye min length: $VIRALFLYE_MIN_LENGTH"
fi
echo ""

# Clean previous results and Nextflow cache
echo "🧹 5. Cleaning previous results and Nextflow cache..."
if [ -d "results_long" ]; then
    echo "Removing previous results directory..."
    rm -rf results_long
fi
# Clean Nextflow work directory to force fresh execution
if [ -d ".nextflow" ]; then
    echo "Clearing Nextflow cache..."
    rm -rf .nextflow/cache
fi

# Set Singularity/Apptainer bind paths and ensure it's enabled
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"
# Ensure Apptainer/Singularity is enabled for Nextflow (don't disable it)
unset NXF_DISABLE_SINGULARITY

# Build nextflow command
NEXTFLOW_CMD="nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db \"$VIRSORTER2_DB\" \
    --deepvirfinder_dir \"$DEEPVIRFINDER_DIR\" \
    --longread true \
    --longread_platform $LONGREAD_PLATFORM"

# Add viralFlye options if enabled
if [ "$ENABLE_VIRALFLYE" = "true" ]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD --enable_viralflye true \
        --viralflye_min_score $VIRALFLYE_MIN_SCORE \
        --viralflye_min_length $VIRALFLYE_MIN_LENGTH"
fi

# Run workflow
echo "🚀 6. Running Metagenome Viral Classification workflow (Long-Read Mode)..."
echo "Command: $NEXTFLOW_CMD"
echo ""
echo "📝 Workflow steps:"
echo "   1. Long-read QC (optional, if enabled)"
echo "   2. metaFlye assembly (--meta mode)"
echo "   3. VirSorter2 viral sequence identification"
echo "   4. DeepVirFinder viral prediction"
echo "   5. Tool result merging (VirSorter2 + DeepVirFinder)"
if [ "$ENABLE_VIRALFLYE" = "true" ]; then
    echo "   6. viralFlye refinement on viral contigs"
    echo "   7. Re-annotation of refined contigs"
    echo "   8. Final merged report for refined contigs"
fi
echo ""
echo "✅ Note: Using metaFlye for long-read assembly"
echo "✅ Note: Platform: $LONGREAD_PLATFORM"
echo "✅ Note: Full run from scratch (no resume mode)"
echo ""

eval $NEXTFLOW_CMD

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $? -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "results_long" ]; then
        echo "📁 Results directory created: results_long/"
        echo "📊 Generated results:"
        
        # Check metaFlye assembly results
        if [ -d "results_long/assembly_metaflye" ]; then
            echo "  ✅ metaFlye assembly: results_long/assembly_metaflye/"
            METAFLYE_CONTIGS=$(find results_long/assembly_metaflye -name "*_metaflye_contigs.fa" | wc -l)
            echo "     - Generated $METAFLYE_CONTIGS contig files"
        fi
        
        # Check VirSorter2 results
        if [ -d "results_long/virsorter2_metaflye" ]; then
            echo "  ✅ VirSorter2 metaFlye results: results_long/virsorter2_metaflye/"
            VS2_METAFLYE=$(find results_long/virsorter2_metaflye -name "*_vs2_final-viral-score.tsv" | wc -l)
            echo "     - Generated $VS2_METAFLYE viral identification reports"
        fi
        
        # Check DeepVirFinder results
        if [ -d "results_long/deepvirfinder_metaflye" ]; then
            echo "  ✅ DeepVirFinder metaFlye results: results_long/deepvirfinder_metaflye/"
            DVF_METAFLYE=$(find results_long/deepvirfinder_metaflye -name "*_dvf_output.txt" | wc -l)
            echo "     - Generated $DVF_METAFLYE viral prediction reports"
        fi
        
        # Check merged viral reports
        if [ -d "results_long/merged_viral_reports_metaflye" ]; then
            echo "  ✅ Merged viral reports (metaFlye): results_long/merged_viral_reports_metaflye/"
            MERGED_VIRAL_METAFLYE=$(find results_long/merged_viral_reports_metaflye -name "*_viral_merged_report.txt" | wc -l)
            echo "     - Generated $MERGED_VIRAL_METAFLYE comprehensive viral analysis reports"
        fi
        
        # Check viralFlye results if enabled
        if [ "$ENABLE_VIRALFLYE" = "true" ]; then
            if [ -d "results_long/assembly_viralflye" ]; then
                echo "  ✅ viralFlye refinement: results_long/assembly_viralflye/"
                VIRALFLYE_CONTIGS=$(find results_long/assembly_viralflye -name "*_viralflye_contigs.fa" | wc -l)
                echo "     - Generated $VIRALFLYE_CONTIGS refined contig files"
            fi
            
            if [ -d "results_long/merged_viral_reports_viralflye" ]; then
                echo "  ✅ Merged viral reports (viralFlye): results_long/merged_viral_reports_viralflye/"
                MERGED_VIRALFLYE=$(find results_long/merged_viral_reports_viralflye -name "*_viralflye_viral_merged_report.txt" | wc -l)
                echo "     - Generated $MERGED_VIRALFLYE refined viral analysis reports"
            fi
        fi
        
        echo ""
        echo "📋 Summary of key viral identification files:"
        echo "  metaFlye assembly contigs:"
        find results_long/assembly_metaflye -name "*_metaflye_contigs.fa" 2>/dev/null | head -5
        echo ""
        echo "  VirSorter2 viral scores:"
        find results_long/virsorter2_metaflye -name "*_vs2_final-viral-score.tsv" 2>/dev/null | head -5
        echo ""
        echo "  DeepVirFinder predictions:"
        find results_long/deepvirfinder_metaflye -name "*_dvf_output.txt" 2>/dev/null | head -5
        echo ""
        echo "  Merged viral reports:"
        find results_long/merged_viral_reports_metaflye -name "*.txt" -o -name "*.csv" 2>/dev/null | head -5
        echo ""
        if [ "$ENABLE_VIRALFLYE" = "true" ]; then
            echo "  ⭐ viralFlye refined contigs:"
            find results_long/assembly_viralflye -name "*_viralflye_contigs.fa" 2>/dev/null | head -5
            echo ""
            echo "  ⭐ viralFlye merged reports:"
            find results_long/merged_viral_reports_viralflye -name "*.txt" -o -name "*.csv" 2>/dev/null | head -5
        fi
        echo ""
        echo "Total files: $(find results_long -type f | wc -l)"
        
    else
        echo "❌ Results directory not found"
    fi
    
else
    echo "❌ Workflow failed with exit code: $?"
    echo "🔍 Check the error log for details"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="

