#!/usr/bin/env nextflow

/*
 * Metagenome Viral Classification Workflow (English Version)
 * 
 * This workflow integrates:
 * 1. Quality control using fastp (optional)
 * 2. Metagenome assembly using MEGAHIT and SPAdes (parallel)
 * 3. Viral sequence identification using VirSorter2 and DeepVirFinder
 * 4. Comprehensive comparison and merging of viral identification results
 * 5. Assembler comparison to identify high-confidence consensus viral sequences
 * 
 * Author: Assistant
 * Version: 5.1.0
 */

nextflow.enable.dsl = 2

// Workflow parameters
// Input data
params.input = null
params.outdir = './results'
params.help = false
// Long-read support
params.longread = false                 // 是否启用长读段分支（PacBio/Nanopore）
params.longread_platform = 'nano'       // 长读段平台: 'nano' 或 'pacbio'
params.skip_longread_qc = true          // 是否跳过长读段质控（默认跳过）
// viralFlye 精修分支
params.enable_viralflye = false         // 是否启用 viralFlye 定向精修（默认关闭）
params.viralflye_min_score = 0.5        // 选择目标 contig 的最小 VirSorter2 分数
params.viralflye_min_length = 1000      // 选择目标 contig 的最小长度（bp）

// Workflow control
params.skip_virsorter2 = false    // Whether to skip VirSorter2 viral identification
params.skip_deepvirfinder = false // Whether to skip DeepVirFinder viral identification
params.skip_merge_reports = false // Whether to skip result merging
params.save_clean_reads = true    // Whether to save filtered clean reads

// Viral identification paths
params.virsorter2_db = null       // VirSorter2 database path
params.deepvirfinder_dir = '/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder' // DeepVirFinder installation directory

// MEGAHIT parameters
params.megahit_memory = 0.8
params.megahit_min_contig_len = 1000

// SPAdes parameters (using metaSPAdes)
params.spades_meta = true

// fastp quality control parameters
params.skip_fastp = false
params.fastp_qualified_quality = 20    // Minimum quality value
params.fastp_unqualified_percent = 40  // Maximum percentage of low-quality bases allowed
params.fastp_min_length = 50           // Minimum read length

// VirSorter2 parameters
params.virsorter2_min_length = 1000      // Minimum contig length
params.virsorter2_min_score = 0.5        // Minimum viral score

// DeepVirFinder parameters  
params.deepvirfinder_min_length = 1000   // Minimum contig length
params.deepvirfinder_pvalue = 0.05       // p-value threshold

// Resource parameters
params.max_cpus = 32
params.max_memory = '512.GB'  // Adjusted to support 512GB memory requirement for SPAdes
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    🦠 Metagenome Viral Classification Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_assembly_classification_workflow_en.nf --input samplesheet.csv --outdir results --virsorter2_db /path/to/db
    
    Required Parameters:
    --input                    Input samplesheet (CSV format)
    --outdir                   Output directory
    --virsorter2_db           VirSorter2 database path
    
    Viral Identification Parameters:
    --deepvirfinder_dir       DeepVirFinder installation directory (default: auto-detected)
    --skip_virsorter2         Skip VirSorter2 analysis (default: false)
    --skip_deepvirfinder      Skip DeepVirFinder analysis (default: false)
    --skip_merge_reports      Skip merging VirSorter2 and DeepVirFinder results (default: false)
    --virsorter2_min_length   Minimum contig length for VirSorter2 (default: 1000)
    --virsorter2_min_score    Minimum viral score for VirSorter2 (default: 0.5)
    --deepvirfinder_min_length Minimum contig length for DeepVirFinder (default: 1000)
    --deepvirfinder_pvalue    P-value threshold for DeepVirFinder (default: 0.05)
    
    Optional Parameters:
    --skip_fastp              Skip fastp quality control (default: false)
    --save_clean_reads        Save filtered clean reads (default: true)
    --longread                Enable long-read mode (PacBio/Nanopore). When true, use metaFlye and skip short-read steps (default: false)
    --longread_platform       Long-read platform: 'nano' (Nanopore) or 'pacbio' (PacBio). Used by Flye (default: nano)
    --skip_longread_qc        Skip long-read QC (default: true)
    --enable_viralflye        Enable viralFlye refinement on viral reads (default: false)
    --viralflye_min_score     Min VS2 score to select targets (default: 0.5)
    --viralflye_min_length    Min contig length to select targets (default: 1000)
    
    Example:
    nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet.csv \\
        --outdir results \\
        --virsorter2_db /scratch/databases/virsorter2/db
    
    Long-read examples:
    - From raw Nanopore/PacBio reads (one FASTQ per sample):
      nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet_long.csv \\
        --outdir results_long \\
        --virsorter2_db /scratch/databases/virsorter2/db \\
        --longread true \\
        --longread_platform nano

    Enable viralFlye refinement:
      ... --enable_viralflye true --viralflye_min_score 0.6 --viralflye_min_length 1500
    
    Long-read samplesheet format (CSV):
      sample,fastq_long
      s1,/path/to/s1_nanopore.fastq.gz
    """
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Input samplesheet is required. Use --input parameter."
}

// Validate VirSorter2 database
if (!params.skip_virsorter2 && !params.virsorter2_db) {
    error "VirSorter2 database path is required. Use --virsorter2_db parameter or --skip_virsorter2 to skip."
}

// Validate DeepVirFinder installation directory
if (!params.skip_deepvirfinder) {
    def dvf_dir = file(params.deepvirfinder_dir)
    if (!dvf_dir.exists() || !dvf_dir.isDirectory()) {
        log.warn "DeepVirFinder directory not found at: ${params.deepvirfinder_dir}"
        log.warn "DeepVirFinder analysis will be skipped."
        params.skip_deepvirfinder = true
    }
}

// Print workflow information
if (!params.longread) {
    log.info """
==========================================
🦠 Metagenome Viral Classification Workflow (Short-Read Mode)
==========================================
Workflow version: 5.1.0
Input samplesheet: ${params.input}
Output directory: ${params.outdir}

Quality Control:
- fastp QC: ${params.skip_fastp ? 'Disabled' : 'Enabled'}
- Save clean reads: ${params.save_clean_reads ? 'Yes' : 'No'}

Assembly Methods:
- MEGAHIT: Enabled
- metaSPAdes: Enabled

Viral Identification:
- VirSorter2: ${params.skip_virsorter2 ? 'Disabled' : 'Enabled'}
${params.skip_virsorter2 ? '' : "  Database: ${params.virsorter2_db}"}
- DeepVirFinder: ${params.skip_deepvirfinder ? 'Disabled' : 'Enabled'}
${params.skip_deepvirfinder ? '' : "  Directory: ${params.deepvirfinder_dir}"}

Result Merging:
- Merge viral reports: ${params.skip_merge_reports ? 'Disabled' : 'Enabled'}
==========================================
"""
} else {
    log.info """
==========================================
🦠 Metagenome Viral Classification Workflow (Long-Read Mode)
==========================================
Workflow version: 5.1.0
Input samplesheet: ${params.input}
Output directory: ${params.outdir}

Long-Read Configuration:
- Platform: ${params.longread_platform} (${params.longread_platform == 'nano' ? 'Nanopore' : 'PacBio'})
- Long-read QC: ${params.skip_longread_qc ? 'Disabled' : 'Enabled'}

Assembly Methods:
- metaFlye: Enabled (--meta mode)

Viral Identification:
- VirSorter2: ${params.skip_virsorter2 ? 'Disabled' : 'Enabled'}
${params.skip_virsorter2 ? '' : "  Database: ${params.virsorter2_db}"}
- DeepVirFinder: ${params.skip_deepvirfinder ? 'Disabled' : 'Enabled'}
${params.skip_deepvirfinder ? '' : "  Directory: ${params.deepvirfinder_dir}"}

viralFlye Refinement:
- Enabled: ${params.enable_viralflye ? 'Yes' : 'No'}
${params.enable_viralflye ? "  Min VS2 score: ${params.viralflye_min_score}" : ''}
${params.enable_viralflye ? "  Min contig length: ${params.viralflye_min_length} bp" : ''}

Result Merging:
- Merge viral reports: ${params.skip_merge_reports ? 'Disabled' : 'Enabled'}
==========================================
"""
}

// 根据模式创建输入通道（短读段或长读段）
if (!params.longread) {
    // 短读段：期望样表包含 fastq_1 与 fastq_2
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def read1 = file(row.fastq_1)
            def read2 = file(row.fastq_2)
            return tuple(sample, [read1, read2])
        }
        .set { ch_reads }
} else {
    // 长读段：期望样表包含 fastq_long（Nanopore 或 PacBio 单端 FASTQ）
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .filter { row ->
            // 过滤空行：检查样本名和路径是否为空
            def sample = row.sample?.trim()
            def read_long = row.fastq_long?.trim()
            return sample && read_long && sample != '' && read_long != ''
        }
        .map { row -> 
            def sample = row.sample?.trim()
            def read_long_path = row.fastq_long?.trim()
            if (!sample || !read_long_path) {
                log.error "Invalid row in samplesheet: sample=${sample}, fastq_long=${read_long_path}"
                return null
            }
            def read_long = file(read_long_path)
            if (!read_long.exists()) {
                log.error "Long-read file not found: ${read_long_path}"
                return null
            }
            return tuple(sample, read_long)
        }
        .filter { it != null }
        .set { ch_long_reads }
    
    // 验证长读段输入通道
    ch_long_reads.view { sample, reads -> "Long-read sample: ${sample}, reads: ${reads}" }
}

// Define workflow
workflow {
    if (!params.longread) {
        // 短读段流程
        // 阶段0：QC（可选）
        if (!params.skip_fastp) {
            FASTP (
                ch_reads
            )
            ch_clean_reads = FASTP.out.clean_reads
        } else {
            ch_clean_reads = ch_reads
        }

        // 阶段1：装配（MEGAHIT + metaSPAdes）
        MEGAHIT_ASSEMBLY (
            ch_clean_reads
        )

        SPADES_ASSEMBLY (
            ch_clean_reads
        )

        // 阶段2：VirSorter2
        if (!params.skip_virsorter2) {
            VIRSORTER2_MEGAHIT (
                MEGAHIT_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )

            VIRSORTER2_SPADES (
                SPADES_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )
        }

        // 阶段3：DeepVirFinder
        if (!params.skip_deepvirfinder) {
            DEEPVIRFINDER_MEGAHIT (
                MEGAHIT_ASSEMBLY.out.contigs
            )

            DEEPVIRFINDER_SPADES (
                SPADES_ASSEMBLY.out.contigs
            )
        }

        // 阶段4：合并与比较（短读段两装配器）
        if (!params.skip_merge_reports && !params.skip_virsorter2 && !params.skip_deepvirfinder) {
            VIRSORTER2_MEGAHIT.out.results
                .join(DEEPVIRFINDER_MEGAHIT.out.results)
                .set { ch_viral_megahit }

            VIRSORTER2_SPADES.out.results
                .join(DEEPVIRFINDER_SPADES.out.results)
                .set { ch_viral_spades }

            MERGE_VIRAL_REPORTS_MEGAHIT (
                ch_viral_megahit
            )

            MERGE_VIRAL_REPORTS_SPADES (
                ch_viral_spades
            )

            // 阶段5：装配器比较
            MERGE_VIRAL_REPORTS_MEGAHIT.out.merged_csv
                .join(MERGE_VIRAL_REPORTS_SPADES.out.merged_csv)
                .set { ch_assembler_comparison }

            COMPARE_ASSEMBLERS (
                ch_assembler_comparison
            )
        }
    } else {
        // 长读段流程
        // 阶段0：长读段QC（可选）
        if (!params.skip_longread_qc) {
            LONGREAD_QC (
                ch_long_reads
            )
            ch_long_clean = LONGREAD_QC.out.clean_long
        } else {
            ch_long_clean = ch_long_reads
        }

        // 阶段1：metaFlye 装配（--meta）
        METAFLYE_ASSEMBLY (
            ch_long_clean
        )

        // 阶段2：VirSorter2（基于 metaFlye 装配的 contigs）
        if (!params.skip_virsorter2) {
            VIRSORTER2_METAFLYE (
                METAFLYE_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )
        }

        // 阶段3：DeepVirFinder（基于 metaFlye 装配的 contigs）
        if (!params.skip_deepvirfinder) {
            DEEPVIRFINDER_METAFLYE (
                METAFLYE_ASSEMBLY.out.contigs
            )
        }

        // 阶段4：合并（单装配器：metaFlye）
        if (!params.skip_merge_reports && !params.skip_virsorter2 && !params.skip_deepvirfinder) {
            VIRSORTER2_METAFLYE.out.results
                .join(DEEPVIRFINDER_METAFLYE.out.results)
                .set { ch_viral_metaflye }

            MERGE_VIRAL_REPORTS_METAFLYE (
                ch_viral_metaflye
            )
        }

        // 阶段5（可选）：viralFlye 精修分支
        if (params.enable_viralflye && !params.skip_virsorter2) {
            // 5.1 选择目标病毒 contigs（基于 VS2 分数与长度阈值）
            SELECT_VIRAL_TARGETS_METAFLYE (
                METAFLYE_ASSEMBLY.out.contigs,
                VIRSORTER2_METAFLYE.out.results
            )

            // 5.2 将长读段与目标 contigs 对齐并抽取相关 reads
            // join 样本名以合并 QC 后的长读段与目标 contigs
            ch_long_clean
                .join(SELECT_VIRAL_TARGETS_METAFLYE.out.targets)
                .set { ch_targets_with_reads }

            SUBSET_LONGREADS_FOR_VIRAL (
                ch_targets_with_reads
            )

            // 5.3 使用（viral）Flye 对选中的 reads 进行定向重装
            VIRALFLYE_REASSEMBLY (
                SUBSET_LONGREADS_FOR_VIRAL.out.selected_reads
            )

            // 5.4 在精修 contigs 上再次进行病毒注释
            if (!params.skip_virsorter2) {
                VIRSORTER2_VIRALFLYE (
                    VIRALFLYE_REASSEMBLY.out.refined_contigs,
                    params.virsorter2_db
                )
            }

            if (!params.skip_deepvirfinder) {
                DEEPVIRFINDER_VIRALFLYE (
                    VIRALFLYE_REASSEMBLY.out.refined_contigs
                )
            }

            if (!params.skip_merge_reports && !params.skip_virsorter2 && !params.skip_deepvirfinder) {
                VIRSORTER2_VIRALFLYE.out.results
                    .join(DEEPVIRFINDER_VIRALFLYE.out.results)
                    .set { ch_viral_viralflye }

                MERGE_VIRAL_REPORTS_VIRALFLYE (
                    ch_viral_viralflye
                )
            }
        }
    }
}

// ================================================================================
// Process Definitions
// ================================================================================

// Process: fastp Quality Control
process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"
    publishDir "${params.outdir}/clean_reads", mode: 'copy', pattern: "*_clean_R*.fastq.gz", enabled: params.save_clean_reads
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_clean_R{1,2}.fastq.gz"), emit: clean_reads
    path("${sample}_fastp.html"), emit: html
    path("${sample}_fastp.json"), emit: json
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    echo "=== fastp Quality Control: ${sample} ==="
    
    # List input files for debugging
    echo "Input files in work directory:"
    ls -lh
    
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_clean_R1.fastq.gz \\
        -O ${sample}_clean_R2.fastq.gz \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
        --length_required ${params.fastp_min_length} \\
        --detect_adapter_for_pe \\
        --compression 6 \\
        --html ${sample}_fastp.html \\
        --json ${sample}_fastp.json
    
    echo "fastp: Quality control completed for ${sample}"
    """
}

// 进程：长读段QC（可选，简化：当前直接拷贝输入，预留对接 Filtlong/NanoFilt）
process LONGREAD_QC {
    tag "${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/longread_qc", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample), path(read_long)

    output:
    tuple val(sample), path("${sample}_long_clean.fastq.gz"), emit: clean_long

    script:
    """
    echo "=== 长读段QC（占位实现）：${sample} ==="
    # 目前默认跳过复杂过滤，仅标准化输出文件名，便于后续流程
    # 如需严格QC，可在此处集成 Filtlong 或 NanoFilt
    if [[ "${read_long}" == *.gz ]]; then
        cp ${read_long} ${sample}_long_clean.fastq.gz
    else
        gzip -c ${read_long} > ${sample}_long_clean.fastq.gz
    fi
    """
}

// Process: MEGAHIT Assembly
process MEGAHIT_ASSEMBLY {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    container 'docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1'
    publishDir "${params.outdir}/assembly_megahit", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_megahit_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== MEGAHIT Assembly: ${sample} ==="
    
    megahit \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o megahit_output \
        -t ${task.cpus} \
        --memory ${params.megahit_memory} \
        --min-contig-len ${params.megahit_min_contig_len}
    
    cp megahit_output/final.contigs.fa ${sample}_megahit_contigs.fa
    
    echo "MEGAHIT: Generated \$(grep -c ">" ${sample}_megahit_contigs.fa) contigs"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    publishDir "${params.outdir}/assembly_spades", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_spades_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== metaSPAdes Assembly: ${sample} ==="
    
    # Use metaSPAdes, disable error correction to avoid memory and bug issues
    metaspades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o spades_output \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()} \
        --only-assembler
    
    cp spades_output/contigs.fasta ${sample}_spades_contigs.fa
    
    echo "metaSPAdes: Generated \$(grep -c ">" ${sample}_spades_contigs.fa) contigs"
    """
}

// 进程：metaFlye 装配（长读段，启用 --meta）
process METAFLYE_ASSEMBLY {
    tag "${sample}_metaFlye"
    label 'process_high'
    conda 'bioconda::flye=2.9'
    publishDir "${params.outdir}/assembly_metaflye", mode: 'copy', pattern: "*.fa"

    input:
    tuple val(sample), path(read_long)

    output:
    tuple val(sample), path("${sample}_metaflye_contigs.fa"), emit: contigs

    script:
    """
    echo "=== metaFlye 装配：${sample} ==="
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_FLAG="--pacbio-raw"
    else
        PLATFORM_FLAG="--nano-raw"
    fi

    flye \\
        \${PLATFORM_FLAG} ${read_long} \\
        --out-dir flye_output \\
        --threads ${task.cpus} \\
        --meta

    cp flye_output/assembly.fasta ${sample}_metaflye_contigs.fa

    echo "metaFlye: 生成 \$(grep -c ">" ${sample}_metaflye_contigs.fa) 条 contigs"
    """
}

// Process: VirSorter2 for MEGAHIT
process VIRSORTER2_MEGAHIT {
    tag "${sample}_MEGAHIT_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_megahit", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_megahit_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    echo "=== VirSorter2 Analysis (MEGAHIT): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_megahit_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_megahit_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_megahit_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_megahit_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from MEGAHIT contigs"
    """
}

// Process: VirSorter2 for SPAdes
process VIRSORTER2_SPADES {
    tag "${sample}_SPAdes_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_spades", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_spades_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_spades_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_spades_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    echo "=== VirSorter2 Analysis (SPAdes): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_spades_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_spades_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_spades_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_spades_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from SPAdes contigs"
    """
}

// 进程：VirSorter2（metaFlye 装配产物）
process VIRSORTER2_METAFLYE {
    tag "${sample}_metaFlye_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'
    publishDir "${params.outdir}/virsorter2_metaflye", mode: 'copy', pattern: "*.{tsv,fa}"

    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)

    output:
    tuple val(sample), path("${sample}_metaflye_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_metaflye_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_metaflye_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true

    script:
    """
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    echo "=== VirSorter2（metaFlye）：${sample} ==="
    virsorter run \
        -i ${contigs} \
        -w virsorter2_output \
        --db-dir ${virsorter2_db} \
        --min-length ${params.virsorter2_min_length} \
        --min-score ${params.virsorter2_min_score} \
        -j ${task.cpus} \
        all

    cp virsorter2_output/final-viral-score.tsv ${sample}_metaflye_vs2_final-viral-score.tsv
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_metaflye_vs2_final-viral-combined.fa
    fi
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_metaflye_vs2_final-viral-boundary.tsv
    fi
    VIRAL_COUNT=\$(tail -n +2 ${sample}_metaflye_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: 识别到 \${VIRAL_COUNT} 条病毒序列（metaFlye）"
    """
}

// Process: DeepVirFinder for MEGAHIT
process DEEPVIRFINDER_MEGAHIT {
    tag "${sample}_MEGAHIT_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_megahit", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (MEGAHIT): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_megahit_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_megahit_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from MEGAHIT contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}

// Process: DeepVirFinder for SPAdes
process DEEPVIRFINDER_SPADES {
    tag "${sample}_SPAdes_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_spades", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (SPAdes): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_spades_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_spades_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from SPAdes contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}

// 进程：DeepVirFinder（metaFlye 装配产物）
process DEEPVIRFINDER_METAFLYE {
    tag "${sample}_metaFlye_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_metaflye", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("${sample}_metaflye_dvf_output.txt"), emit: results

    script:
    """
    echo "=== DeepVirFinder（metaFlye）：${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    set +u
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf 环境未找到: \$DVF_ENV"; exit 1; fi
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    [ -z "\$CONDA_BASE" ] && CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then source "\$CONDA_BASE/etc/profile.d/conda.sh"; else echo "❌ 无法找到 conda.sh"; exit 1; fi
    conda activate "\$DVF_ENV" || { echo "❌ 激活 dvf 失败"; exit 1; }
    export PATH="\$DVF_ENV/bin:\$PATH"; export CONDA_PREFIX="\$DVF_ENV"; export CONDA_DEFAULT_ENV="dvf"
    unset PYTHONPATH
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    export KERAS_BACKEND=theano
    set -u
    python ${params.deepvirfinder_dir}/dvf.py \
        -i ${contigs} \
        -o dvf_output \
        -l ${params.deepvirfinder_min_length} \
        -c ${task.cpus}
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_metaflye_dvf_output.txt
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_metaflye_dvf_output.txt)
    echo "DeepVirFinder: 预测到 \${VIRAL_COUNT} 条病毒序列（metaFlye，p<${params.deepvirfinder_pvalue}})"
    """
}
// Process: Merge Viral Identification Reports for MEGAHIT
// Integrate VirSorter2 and DeepVirFinder results
process MERGE_VIRAL_REPORTS_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_megahit_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_megahit_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_megahit_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_virsorter2(file_path):
        \"\"\"
        Parse VirSorter2 output file
        Columns: seqname, dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae, max_score, max_score_group, length, hallmark, viral_gene, cellular_gene
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            # Extract sequence names and scores
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['seqname']
                # Normalize sequence name: remove ||full, ||partial, etc.
                seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
                max_score = row['max_score']
                viral_dict[seqname_normalized] = {
                    'vs2_score': max_score,
                    'vs2_group': row['max_score_group'],
                    'vs2_length': row['length']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse VirSorter2 results: {e}")
            return {}
    
    def parse_deepvirfinder(file_path):
        \"\"\"
        Parse DeepVirFinder output file
        Columns: name, len, score, pvalue
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['name']
                # Normalize sequence name: remove flags and metadata (keep only the contig ID)
                seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
                viral_dict[seqname_normalized] = {
                    'dvf_score': row['score'],
                    'dvf_pvalue': row['pvalue'],
                    'dvf_length': row['len']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse DeepVirFinder results: {e}")
            return {}
    
    # Parse result files
    print(f"Parsing VirSorter2 results: ${virsorter2_results}")
    vs2_dict = parse_virsorter2("${virsorter2_results}")
    print(f"VirSorter2: Parsed {len(vs2_dict)} sequences")
    if len(vs2_dict) > 0:
        print(f"Sample VirSorter2 sequences: {list(vs2_dict.keys())[:5]}")
    
    print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
    dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
    print(f"DeepVirFinder: Parsed {len(dvf_dict)} sequences")
    if len(dvf_dict) > 0:
        print(f"Sample DeepVirFinder sequences: {list(dvf_dict.keys())[:5]}")
    
    # Merge results
    all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
    print(f"Total unique sequences: {len(all_sequences)}")
    if len(all_sequences) == 0:
        print("⚠️ WARNING: No viral sequences found at all! Please check:")
        print("  1. Are there any sequences in the VirSorter2 output file?")
        print("  2. Are there any sequences in the DeepVirFinder output file?")
        print("  3. Check the detection thresholds in the config file")
    
    # Generate comprehensive report
    with open("${sample}_megahit_viral_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Viral Identification Comprehensive Analysis Report - MEGAHIT Assembly Results\\n")
        f.write("VirSorter2 + DeepVirFinder\\n")
        f.write("="*80 + "\\n\\n")
        
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
        f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    
        # Count consensus sequences (identified by both tools)
        consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
        f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    
        # Identified by only one tool
        vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
        dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
        f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
        f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")
    
        # DVF significant sequences (p-value < ${params.deepvirfinder_pvalue})
        dvf_significant = [seq for seq, data in dvf_dict.items() 
                          if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
        f.write(f"DeepVirFinder significant sequences (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
    
        # Consensus sequence details (recommended high-confidence viral sequences)
        f.write("\\n[Consensus Viral Sequences (High Confidence)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for seq in sorted(consensus):
            vs2_score = vs2_dict[seq]['vs2_score']
            dvf_score = dvf_dict[seq]['dvf_score']
            dvf_pval = dvf_dict[seq]['dvf_pvalue']
            f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save detailed data in CSV format
    merged_data = []
    for seq in all_sequences:
        row = {
            'sequence_name': seq,
            'identified_by': '',
            'vs2_score': None,
            'vs2_group': None,
            'dvf_score': None,
            'dvf_pvalue': None,
            'consensus': False
        }
        
        if seq in vs2_dict:
            row['vs2_score'] = vs2_dict[seq]['vs2_score']
            row['vs2_group'] = vs2_dict[seq]['vs2_group']
        
        if seq in dvf_dict:
            row['dvf_score'] = dvf_dict[seq]['dvf_score']
            row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
        
        if seq in consensus:
            row['identified_by'] = 'Both'
            row['consensus'] = True
        elif seq in vs2_only:
            row['identified_by'] = 'VirSorter2_only'
        else:
            row['identified_by'] = 'DeepVirFinder_only'
        
        merged_data.append(row)
    
    merged_df = pd.DataFrame(merged_data)
    merged_df.to_csv("${sample}_megahit_viral_merged_report.csv", index=False)
    
    # Save consensus sequence list (recommended for downstream analysis)
    with open("${sample}_megahit_viral_consensus.txt", 'w') as f:
        for seq in sorted(consensus):
            f.write(seq + "\\n")
    
    print(f"Viral identification report generated successfully: ${sample} (MEGAHIT)")
    print(f"Consensus viral sequences: {len(consensus)}")
    """
}

// Process: Merge Viral Identification Reports for SPAdes
process MERGE_VIRAL_REPORTS_SPADES {
    tag "${sample}_SPAdes"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_spades_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_spades_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_spades_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_virsorter2(file_path):
        \"\"\"
        Parse VirSorter2 output file
        Columns: seqname, dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae, max_score, max_score_group, length, hallmark, viral_gene, cellular_gene
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            # Extract sequence names and scores
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['seqname']
                # Normalize sequence name: remove ||full, ||partial, etc.
                seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
                max_score = row['max_score']
                viral_dict[seqname_normalized] = {
                    'vs2_score': max_score,
                    'vs2_group': row['max_score_group'],
                    'vs2_length': row['length']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse VirSorter2 results: {e}")
            return {}
    
    def parse_deepvirfinder(file_path):
        \"\"\"
        Parse DeepVirFinder output file
        Columns: name, len, score, pvalue
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['name']
                # Normalize sequence name: remove flags and metadata (keep only the contig ID)
                seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
                viral_dict[seqname_normalized] = {
                    'dvf_score': row['score'],
                    'dvf_pvalue': row['pvalue'],
                    'dvf_length': row['len']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse DeepVirFinder results: {e}")
            return {}
    
    # Parse result files
    print(f"Parsing VirSorter2 results: ${virsorter2_results}")
    vs2_dict = parse_virsorter2("${virsorter2_results}")
    print(f"VirSorter2: Parsed {len(vs2_dict)} sequences")
    if len(vs2_dict) > 0:
        print(f"Sample VirSorter2 sequences: {list(vs2_dict.keys())[:5]}")
    
    print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
    dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
    print(f"DeepVirFinder: Parsed {len(dvf_dict)} sequences")
    if len(dvf_dict) > 0:
        print(f"Sample DeepVirFinder sequences: {list(dvf_dict.keys())[:5]}")
    
    # Merge results
    all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
    print(f"Total unique sequences: {len(all_sequences)}")
    if len(all_sequences) == 0:
        print("⚠️ WARNING: No viral sequences found at all! Please check:")
        print("  1. Are there any sequences in the VirSorter2 output file?")
        print("  2. Are there any sequences in the DeepVirFinder output file?")
        print("  3. Check the detection thresholds in the config file")
    
    # Generate comprehensive report
    with open("${sample}_spades_viral_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Viral Identification Comprehensive Analysis Report - SPAdes Assembly Results\\n")
        f.write("VirSorter2 + DeepVirFinder\\n")
        f.write("="*80 + "\\n\\n")
    
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
        f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    
        # Count consensus sequences (identified by both tools)
        consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
        f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    
        # Identified by only one tool
        vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
        dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
        f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
        f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")
    
        # DVF significant sequences (p-value < ${params.deepvirfinder_pvalue})
        dvf_significant = [seq for seq, data in dvf_dict.items() 
                          if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
        f.write(f"DeepVirFinder significant sequences (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
    
        # Consensus sequence details (recommended high-confidence viral sequences)
        f.write("\\n[Consensus Viral Sequences (High Confidence)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for seq in sorted(consensus):
            vs2_score = vs2_dict[seq]['vs2_score']
            dvf_score = dvf_dict[seq]['dvf_score']
            dvf_pval = dvf_dict[seq]['dvf_pvalue']
            f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save detailed data in CSV format
    merged_data = []
    for seq in all_sequences:
        row = {
            'sequence_name': seq,
            'identified_by': '',
            'vs2_score': None,
            'vs2_group': None,
            'dvf_score': None,
            'dvf_pvalue': None,
            'consensus': False
        }
        
        if seq in vs2_dict:
            row['vs2_score'] = vs2_dict[seq]['vs2_score']
            row['vs2_group'] = vs2_dict[seq]['vs2_group']
        
        if seq in dvf_dict:
            row['dvf_score'] = dvf_dict[seq]['dvf_score']
            row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
        
        if seq in consensus:
            row['identified_by'] = 'Both'
            row['consensus'] = True
        elif seq in vs2_only:
            row['identified_by'] = 'VirSorter2_only'
        else:
            row['identified_by'] = 'DeepVirFinder_only'
        
        merged_data.append(row)
    
    merged_df = pd.DataFrame(merged_data)
    merged_df.to_csv("${sample}_spades_viral_merged_report.csv", index=False)
    
    # Save consensus sequence list (recommended for downstream analysis)
    with open("${sample}_spades_viral_consensus.txt", 'w') as f:
        for seq in sorted(consensus):
            f.write(seq + "\\n")
    
    print(f"Viral identification report generated successfully: ${sample} (SPAdes)")
    print(f"Consensus viral sequences: {len(consensus)}")
    """
}

// 进程：合并（metaFlye）——整合 VirSorter2 与 DeepVirFinder 结果
process MERGE_VIRAL_REPORTS_METAFLYE {
    tag "${sample}_metaFlye"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_metaflye", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)

    output:
    tuple val(sample), path("${sample}_metaflye_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_metaflye_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_metaflye_viral_consensus.txt"), emit: consensus_list

    script:
    """
    python3 << 'PYTHON_SCRIPT'
import pandas as pd

def parse_virsorter2(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname']
            seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'vs2_score': row['max_score'],
                'vs2_group': row['max_score_group'],
                'vs2_length': row['length']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse VirSorter2 results: {e}")
        return {}

def parse_deepvirfinder(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['name']
            seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'dvf_score': row['score'],
                'dvf_pvalue': row['pvalue'],
                'dvf_length': row['len']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse DeepVirFinder results: {e}")
        return {}

vs2_dict = parse_virsorter2("${virsorter2_results}")
dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
with open("${sample}_metaflye_viral_merged_report.txt", 'w', encoding='utf-8') as f:
    f.write("="*80 + "\\n")
    f.write("Viral Identification Comprehensive Analysis Report - metaFlye Assembly Results\\n")
    f.write("VirSorter2 + DeepVirFinder\\n")
    f.write("="*80 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
    f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
    f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")

merged_data = []
for seq in all_sequences:
    row = {
        'sequence_name': seq,
        'identified_by': '',
        'vs2_score': None,
        'vs2_group': None,
        'dvf_score': None,
        'dvf_pvalue': None,
        'consensus': False
    }
    if seq in vs2_dict:
        row['vs2_score'] = vs2_dict[seq]['vs2_score']
        row['vs2_group'] = vs2_dict[seq]['vs2_group']
    if seq in dvf_dict:
        row['dvf_score'] = dvf_dict[seq]['dvf_score']
        row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
    if seq in consensus:
        row['identified_by'] = 'Both'
        row['consensus'] = True
    elif seq in vs2_only:
        row['identified_by'] = 'VirSorter2_only'
    else:
        row['identified_by'] = 'DeepVirFinder_only'
    merged_data.append(row)

pd.DataFrame(merged_data).to_csv("${sample}_metaflye_viral_merged_report.csv", index=False)

with open("${sample}_metaflye_viral_consensus.txt", 'w') as f:
    for seq in sorted(consensus):
        f.write(seq + "\\n")

print(f"Viral identification report generated successfully: ${sample} (metaFlye)")
PYTHON_SCRIPT
    """
}

// 进程：选择 VS2 阈值以上的目标病毒 contigs，并导出子集 FASTA
process SELECT_VIRAL_TARGETS_METAFLYE {
    tag "${sample}_select_targets"
    label 'process_low'
    conda 'bioconda::seqkit=2.8.2 pandas'
    publishDir "${params.outdir}/viralflye_targets", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(contigs)
    tuple val(sample2), path(vs2_scores)

    output:
    tuple val(sample), path("${sample}_viral_targets.fa"), emit: targets
    path("${sample}_viral_target_ids.txt"), emit: ids

    when:
    sample == sample2

    script:
    """
    echo "=== 选择病毒目标 contigs：${sample} ==="
    python - << 'PY'
import pandas as pd, sys
vs2 = pd.read_csv("${vs2_scores}", sep='\t')
ok = (vs2['max_score'] >= ${params.viralflye_min_score}) & (vs2['length'] >= ${params.viralflye_min_length})
ids = []
for _, r in vs2[ok].iterrows():
    name = r['seqname']
    name = name.split('||')[0] if '||' in name else name
    ids.append(name)
open("${sample}_viral_target_ids.txt", 'w').write('\n'.join(ids)+'\n')
PY
    # 用 seqkit 提取目标 contigs
    if [ -s ${sample}_viral_target_ids.txt ]; then
        seqkit grep -f ${sample}_viral_target_ids.txt ${contigs} > ${sample}_viral_targets.fa
    else
        # 若无目标，生成空文件避免中断
        echo > ${sample}_viral_targets.fa
    fi
    """
}

// 进程：对目标 contigs 进行 minimap2 比对并用 samtools 抽取相关 reads
process SUBSET_LONGREADS_FOR_VIRAL {
    tag "${sample}_subset_reads"
    label 'process_medium'
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'
    publishDir "${params.outdir}/viralflye_reads", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample), path(read_long), path(targets_fa)

    output:
    tuple val(sample), path("${sample}_viral_reads.fastq.gz"), emit: selected_reads

    script:
    """
    echo "=== 抽取病毒相关 reads：${sample} ==="
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_OPT="map-pb"
    else
        PLATFORM_OPT="map-ont"
    fi
    # 生成 SAM 并提取比对上的 reads 到 FASTQ
    minimap2 -t ${task.cpus} -x \${PLATFORM_OPT} -a ${targets_fa} ${read_long} \\
      | samtools view -b -F 4 -@ ${task.cpus} \\
      | samtools fastq -@ ${task.cpus} -n - \\
      | gzip -c > ${sample}_viral_reads.fastq.gz
    """
}

// 进程：viralFlye 定向重装（此处使用 Flye 实现定向 reassembly）
process VIRALFLYE_REASSEMBLY {
    tag "${sample}_viralFlye"
    label 'process_high'
    conda 'bioconda::flye=2.9'
    publishDir "${params.outdir}/assembly_viralflye", mode: 'copy', pattern: "*.fa"

    input:
    tuple val(sample), path(viral_reads)

    output:
    tuple val(sample), path("${sample}_viralflye_contigs.fa"), emit: refined_contigs

    script:
    """
    echo "=== viralFlye 定向重装：${sample} ==="
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_FLAG="--pacbio-raw"
    else
        PLATFORM_FLAG="--nano-raw"
    fi
    flye \${PLATFORM_FLAG} ${viral_reads} --out-dir viralflye_out --threads ${task.cpus} --meta
    cp viralflye_out/assembly.fasta ${sample}_viralflye_contigs.fa
    echo "viralFlye: 生成 \$(grep -c ">" ${sample}_viralflye_contigs.fa) 条 contigs"
    """
}

// 进程：VirSorter2（viralFlye 精修 contigs）
process VIRSORTER2_VIRALFLYE {
    tag "${sample}_viralFlye_VS2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'
    publishDir "${params.outdir}/virsorter2_viralflye", mode: 'copy', pattern: "*.{tsv,fa}"

    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)

    output:
    tuple val(sample), path("${sample}_viralflye_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_viralflye_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_viralflye_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true

    script:
    """
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    echo "=== VirSorter2（viralFlye）：${sample} ==="
    virsorter run -i ${contigs} -w vs2_out --db-dir ${virsorter2_db} \
        --min-length ${params.virsorter2_min_length} --min-score ${params.virsorter2_min_score} -j ${task.cpus} all
    cp vs2_out/final-viral-score.tsv ${sample}_viralflye_vs2_final-viral-score.tsv
    if [ -f vs2_out/final-viral-combined.fa ]; then cp vs2_out/final-viral-combined.fa ${sample}_viralflye_vs2_final-viral-combined.fa; fi
    if [ -f vs2_out/final-viral-boundary.tsv ]; then cp vs2_out/final-viral-boundary.tsv ${sample}_viralflye_vs2_final-viral-boundary.tsv; fi
    """
}

// 进程：DeepVirFinder（viralFlye 精修 contigs）
process DEEPVIRFINDER_VIRALFLYE {
    tag "${sample}_viralFlye_DVF"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_viralflye", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("${sample}_viralflye_dvf_output.txt"), emit: results

    script:
    """
    echo "=== DeepVirFinder（viralFlye）：${sample} ==="
    set +u
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    [ -d "\$DVF_ENV" ] || { echo "❌ dvf 环境未找到: \$DVF_ENV"; exit 1; }
    CONDA_BASE=\$(conda info --base 2>/dev/null); [ -z "\$CONDA_BASE" ] && CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ] && source "\$CONDA_BASE/etc/profile.d/conda.sh" || { echo "❌ 无法找到 conda.sh"; exit 1; }
    conda activate "\$DVF_ENV" || { echo "❌ 激活 dvf 失败"; exit 1; }
    export PATH="\$DVF_ENV/bin:\$PATH"; unset PYTHONPATH; export KERAS_BACKEND=theano
    set -u
    python ${params.deepvirfinder_dir}/dvf.py -i ${contigs} -o dvf_out -l ${params.deepvirfinder_min_length} -c ${task.cpus}
    cp dvf_out/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_viralflye_dvf_output.txt
    """
}

// 进程：合并（viralFlye 精修产物）
process MERGE_VIRAL_REPORTS_VIRALFLYE {
    tag "${sample}_viralFlye"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_viralflye", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)

    output:
    tuple val(sample), path("${sample}_viralflye_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_viralflye_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_viralflye_viral_consensus.txt"), emit: consensus_list

    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    import pandas as pd
    def parse_vs2(p):
        try:
            df = pd.read_csv(p, sep='\t'); d={}
            for _,r in df.iterrows():
                n=r['seqname'].split('||')[0] if '||' in r['seqname'] else r['seqname']
                d[n]={'vs2_score':r['max_score'],'vs2_group':r['max_score_group'],'vs2_length':r['length']}
            return d
        except Exception as e:
            print('Warn VS2',e); return {}
    def parse_dvf(p):
        try:
            df = pd.read_csv(p, sep='\t'); d={}
            for _,r in df.iterrows():
                n=r['name'].split()[0] if ' ' in r['name'] else r['name']
                d[n]={'dvf_score':r['score'],'dvf_pvalue':r['pvalue'],'dvf_length':r['len']}
            return d
        except Exception as e:
            print('Warn DVF',e); return {}
    vs2=parse_vs2("${virsorter2_results}"); dvf=parse_dvf("${deepvirfinder_results}")
    all_ids=set(vs2)|set(dvf); consensus=set(vs2)&set(dvf)
    rows=[]
    for s in all_ids:
        row={'sequence_name':s,'identified_by':'','vs2_score':None,'vs2_group':None,'dvf_score':None,'dvf_pvalue':None,'consensus':False}
        if s in vs2: row['vs2_score']=vs2[s]['vs2_score']; row['vs2_group']=vs2[s]['vs2_group']
        if s in dvf: row['dvf_score']=dvf[s]['dvf_score']; row['dvf_pvalue']=dvf[s]['dvf_pvalue']
        if s in consensus: row['identified_by']='Both'; row['consensus']=True
        elif s in set(vs2)-set(dvf): row['identified_by']='VirSorter2_only'
        else: row['identified_by']='DeepVirFinder_only'
        rows.append(row)
    pd.DataFrame(rows).to_csv("${sample}_viralflye_viral_merged_report.csv", index=False)
    with open("${sample}_viralflye_viral_merged_report.txt",'w') as f:
        f.write('Viral report (viralFlye)\n')
    with open("${sample}_viralflye_viral_consensus.txt",'w') as f:
        for s in sorted(consensus): f.write(s+'\n')
    print(f"viralFlye merged report written: ${sample}")
    """
}

// ================================================================================
// Assembler Comparison
// ================================================================================

// Process: Compare MEGAHIT vs SPAdes Viral Identification Results
process COMPARE_ASSEMBLERS {
    tag "${sample}_Assembler_Comparison"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/assembler_comparison", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(megahit_report), path(spades_report)
    
    output:
    tuple val(sample), path("${sample}_assembler_comparison.txt"), emit: comparison_report
    path("${sample}_assembler_comparison.csv"), emit: comparison_csv
    path("${sample}_consensus_viral_sequences.txt"), emit: final_consensus
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_merged_report(file_path):
        \"\"\"Parse merged viral report CSV file\"\"\"
        import os
        try:
            # If input is already CSV, use it directly
            # Otherwise try to convert .txt to .csv
            if file_path.endswith('.csv'):
                csv_path = file_path
            else:
                csv_path = file_path.replace('.txt', '.csv')
            
            print(f"  Attempting to read CSV from: {csv_path}")
            if not os.path.exists(csv_path):
                print(f"  ⚠️  CSV file does not exist: {csv_path}")
                print(f"  Checking if original file exists: {file_path}")
                if os.path.exists(file_path):
                    print(f"  Original file exists, but CSV version not found")
                return pd.DataFrame()
            
            df = pd.read_csv(csv_path)
            print(f"  ✅ Successfully parsed {len(df)} rows from {csv_path}")
            return df
        except Exception as e:
            print(f"  ❌ Warning: Failed to parse {file_path}: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    print(f"\\nParsing MEGAHIT results from: ${megahit_report}")
    print(f"  Input file path: ${megahit_report}")
    
    megahit_df = parse_merged_report("${megahit_report}")
    print(f"MEGAHIT DataFrame shape: {megahit_df.shape}")
    print(f"MEGAHIT columns: {list(megahit_df.columns)}")
    if len(megahit_df) > 0:
        print(f"MEGAHIT first 3 rows:\\n{megahit_df.head(3)}")
    else:
        print("⚠️  WARNING: MEGAHIT DataFrame is empty!")
    
    print(f"\\nParsing SPAdes results from: ${spades_report}")
    print(f"  Input file path: ${spades_report}")
    
    spades_df = parse_merged_report("${spades_report}")
    print(f"SPAdes DataFrame shape: {spades_df.shape}")
    print(f"SPAdes columns: {list(spades_df.columns)}")
    if len(spades_df) > 0:
        print(f"SPAdes first 3 rows:\\n{spades_df.head(3)}")
    else:
        print("⚠️  WARNING: SPAdes DataFrame is empty!")
    
    # Extract viral sequences (only "Both" type - high confidence sequences)
    print(f"\\nExtracting 'Both' sequences from MEGAHIT...")
    print(f"MEGAHIT columns: {list(megahit_df.columns) if len(megahit_df) > 0 else 'empty'}")
    if len(megahit_df) > 0 and 'identified_by' in megahit_df.columns:
        print(f"MEGAHIT 'identified_by' value counts:\\n{megahit_df['identified_by'].value_counts()}")
        megahit_both_df = megahit_df[megahit_df['identified_by'] == 'Both']
        print(f"MEGAHIT 'Both' sequences found: {len(megahit_both_df)}")
        megahit_seqs = set(megahit_both_df['sequence_name'].tolist())
        print(f"MEGAHIT sequences extracted: {len(megahit_seqs)}")
    else:
        print(f"⚠️  MEGAHIT extraction failed - columns: {list(megahit_df.columns) if len(megahit_df) > 0 else 'empty'}")
        megahit_seqs = set()
    
    print(f"\\nExtracting 'Both' sequences from SPAdes...")
    print(f"SPAdes columns: {list(spades_df.columns) if len(spades_df) > 0 else 'empty'}")
    if len(spades_df) > 0 and 'identified_by' in spades_df.columns:
        print(f"SPAdes 'identified_by' value counts:\\n{spades_df['identified_by'].value_counts()}")
        spades_both_df = spades_df[spades_df['identified_by'] == 'Both']
        print(f"SPAdes 'Both' sequences found: {len(spades_both_df)}")
        spades_seqs = set(spades_both_df['sequence_name'].tolist())
        print(f"SPAdes sequences extracted: {len(spades_seqs)}")
    else:
        print(f"⚠️  SPAdes extraction failed - columns: {list(spades_df.columns) if len(spades_df) > 0 else 'empty'}")
        spades_seqs = set()
    
    print(f"\\nMEGAHIT 'Both' sequences (high confidence): {len(megahit_seqs)}")
    print(f"SPAdes 'Both' sequences (high confidence): {len(spades_seqs)}")
    
    # Note: Sequence IDs from different assemblers are different (k141_XXX vs NODE_XXX)
    # "Consensus" means both assemblers detected viruses, not identical sequence IDs
    # We output all high-confidence sequences from both assemblers
    megahit_only = megahit_seqs
    spades_only = spades_seqs
    all_viral = megahit_seqs | spades_seqs
    
    # Consensus: if both assemblers detected viruses, all sequences are "consensus"
    # Otherwise, consensus is empty
    if len(megahit_seqs) > 0 and len(spades_seqs) > 0:
        consensus_both_assemblers = all_viral
        print(f"Both assemblers detected viruses: {len(megahit_seqs)} (MEGAHIT) + {len(spades_seqs)} (SPAdes) = {len(all_viral)} total consensus sequences")
    else:
        consensus_both_assemblers = set()
        print(f"Only one assembler detected viruses: {len(megahit_seqs)} (MEGAHIT) + {len(spades_seqs)} (SPAdes)")
    consensus_count = len(consensus_both_assemblers)
    
    # Generate comprehensive comparison report
    with open("${sample}_assembler_comparison.txt", 'w', encoding='utf-8') as f:
        f.write("="*100 + "\\n")
        f.write("Assembler Comparison Report - Viral Identification Results\\n")
        f.write("MEGAHIT vs metaSPAdes\\n")
        f.write("Sample: ${sample}\\n")
        f.write("="*100 + "\\n\\n")
        
        f.write("[Overall Statistics]\\n")
        f.write("-"*100 + "\\n")
        f.write(f"MEGAHIT identified viral sequences:    {len(megahit_seqs):,}\\n")
        f.write(f"SPAdes identified viral sequences:     {len(spades_seqs):,}\\n")
        f.write(f"Total unique viral sequences:          {len(all_viral):,}\\n\\n")
        
        f.write(f"Consensus viral sequences (both assemblies detected): {len(consensus_both_assemblers):,}\\n")
        f.write(f"MEGAHIT viral sequences:                {len(megahit_seqs):,}\\n")
        f.write(f"SPAdes viral sequences:                {len(spades_seqs):,}\\n\\n")
        
        # Calculate consistency
        if len(all_viral) > 0:
            consistency = len(consensus_both_assemblers) / len(all_viral) * 100
            f.write(f"Assembler consistency:                 {consistency:.2f}%\\n\\n")
        
        f.write("="*100 + "\\n")
        f.write("[Recommendation]\\n")
        f.write("-"*100 + "\\n")
        f.write(f"High-confidence viral sequences (identified by both methods): {len(consensus_both_assemblers):,}\\n")
        f.write("Recommend prioritizing these consensus sequences for downstream analysis.\\n\\n")
        
        f.write("[Detailed Analysis]\\n")
        f.write("-"*100 + "\\n")
        
        # MEGAHIT advantages
        if len(megahit_only) > 0:
            f.write(f"\\nMEGAHIT-specific sequences ({len(megahit_only)}):\\n")
            f.write("  - May represent low-coverage or high-complexity regions\\n")
            f.write("  - MEGAHIT has stronger assembly capability for complex structures\\n")
        
        # SPAdes advantages
        if len(spades_only) > 0:
            f.write(f"\\nSPAdes-specific sequences ({len(spades_only)}):\\n")
            f.write("  - May represent high-coverage regions\\n")
            f.write("  - SPAdes kmer strategy may capture more details\\n")
        
        f.write("\\n" + "="*100 + "\\n")
        f.write("[Statistical Summary]\\n")
        f.write("-"*100 + "\\n")
        
        # Calculate consensus ratio for each assembler
        if len(megahit_seqs) > 0:
            megahit_consensus_pct = len(consensus_both_assemblers) / len(megahit_seqs) * 100
            f.write(f"Consensus ratio in MEGAHIT sequences:  {megahit_consensus_pct:.2f}%\\n")
        
        if len(spades_seqs) > 0:
            spades_consensus_pct = len(consensus_both_assemblers) / len(spades_seqs) * 100
            f.write(f"Consensus ratio in SPAdes sequences:   {spades_consensus_pct:.2f}%\\n")
    
    # Generate CSV detailed comparison
    # Note: Since sequence IDs are different between assemblers, we treat all as consensus if both assemblers detected viruses
    comparison_data = []
    
    for seq in all_viral:
        row = {
            'sequence': seq,
            'found_in_MEGAHIT': 'Yes' if seq in megahit_seqs else 'No',
            'found_in_SPAdes': 'Yes' if seq in spades_seqs else 'No',
            'status': 'Consensus' if len(consensus_both_assemblers) > 0 else 'Single_assembler'
        }
        
        # Add MEGAHIT information
        if seq in megahit_seqs:
            megahit_row = megahit_df[megahit_df['sequence_name'] == seq].iloc[0]
            row['MEGAHIT_vs2_score'] = megahit_row.get('vs2_score', 'N/A')
            row['MEGAHIT_dvf_score'] = megahit_row.get('dvf_score', 'N/A')
            row['MEGAHIT_identified_by'] = megahit_row.get('identified_by', 'N/A')
        else:
            row['MEGAHIT_vs2_score'] = 'N/A'
            row['MEGAHIT_dvf_score'] = 'N/A'
            row['MEGAHIT_identified_by'] = 'N/A'
        
        # Add SPAdes information
        if seq in spades_seqs:
            spades_row = spades_df[spades_df['sequence_name'] == seq].iloc[0]
            row['SPAdes_vs2_score'] = spades_row.get('vs2_score', 'N/A')
            row['SPAdes_dvf_score'] = spades_row.get('dvf_score', 'N/A')
            row['SPAdes_identified_by'] = spades_row.get('identified_by', 'N/A')
        else:
            row['SPAdes_vs2_score'] = 'N/A'
            row['SPAdes_dvf_score'] = 'N/A'
            row['SPAdes_identified_by'] = 'N/A'
        
        comparison_data.append(row)
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df.to_csv("${sample}_assembler_comparison.csv", index=False)
    
    # Save final consensus sequence list
    with open("${sample}_consensus_viral_sequences.txt", 'w') as f:
        if len(consensus_both_assemblers) > 0:
            f.write("# High-confidence viral sequences from both MEGAHIT and SPAdes\\n")
            f.write("# Both assemblers detected viruses (consensus by existence)\\n")
        else:
            f.write("# Viral sequences from single assembler only\\n")
        f.write(f"# Sample: ${sample}\\n")
        f.write(f"# Total sequences: {len(consensus_both_assemblers)}\\n")
        f.write("# Note: Sequence IDs differ between assemblers, consensus means both detected viruses\\n")
        f.write("#\\n")
        for seq in sorted(consensus_both_assemblers):
            f.write(seq + "\\n")
    
    print(f"\\nAssembler comparison complete: ${sample}")
    print(f"  MEGAHIT: {len(megahit_seqs)} viral sequences")
    print(f"  SPAdes:  {len(spades_seqs)} viral sequences")
    print(f"  Consensus: {len(consensus_both_assemblers)} viral sequences")
    if len(consensus_both_assemblers) > 0:
        print(f"  Both assemblers detected viruses: {len(megahit_seqs)} + {len(spades_seqs)} = {len(consensus_both_assemblers)} total consensus sequences")
    """
}

// Workflow completion message
workflow.onComplete {
    log.info """
    ==========================================
    🎯 Metagenome Viral Classification Results
    ==========================================
    Pipeline completed successfully!
    
    Results directory: ${params.outdir}
    
    Generated files:
    - fastp/: Quality control reports
      * *_fastp.html: HTML quality reports
      * *_fastp.json: JSON quality data
      
    - clean_reads/: Filtered clean reads (if save_clean_reads=true)
      * *_clean_R1.fastq.gz: Forward clean reads
      * *_clean_R2.fastq.gz: Reverse clean reads
      
    - assembly_megahit/: MEGAHIT assembly results
      * *_megahit_contigs.fa: Assembled contigs
      
    - assembly_spades/: metaSPAdes assembly results
      * *_spades_contigs.fa: Assembled contigs
      
    - virsorter2_megahit/: VirSorter2 viral identification (MEGAHIT)
      * *_megahit_vs2_final-viral-score.tsv: Viral scores
      * *_megahit_vs2_final-viral-combined.fa: Identified viral contigs
      
    - virsorter2_spades/: VirSorter2 viral identification (SPAdes)
      * *_spades_vs2_final-viral-score.tsv: Viral scores
      * *_spades_vs2_final-viral-combined.fa: Identified viral contigs
      
    - deepvirfinder_megahit/: DeepVirFinder viral prediction (MEGAHIT)
      * *_megahit_dvf_output.txt: Prediction results with scores and p-values
      
    - deepvirfinder_spades/: DeepVirFinder viral prediction (SPAdes)
      * *_spades_dvf_output.txt: Prediction results with scores and p-values
      
    - merged_viral_reports_megahit/: Integrated viral analysis (MEGAHIT)
      * *_megahit_viral_merged_report.txt: Comprehensive viral identification report
      * *_megahit_viral_merged_report.csv: Detailed comparison data
      * *_megahit_viral_consensus.txt: High-confidence viral sequences list
      
    - merged_viral_reports_spades/: Integrated viral analysis (SPAdes)
      * *_spades_viral_merged_report.txt: Comprehensive viral identification report
      * *_spades_viral_merged_report.csv: Detailed comparison data
      * *_spades_viral_consensus.txt: High-confidence viral sequences list
    
    - assembler_comparison/: MEGAHIT vs SPAdes comparison ⭐
      * *_assembler_comparison.txt: Comprehensive assembler comparison report
      * *_assembler_comparison.csv: Detailed comparison data
      * *_consensus_viral_sequences.txt: Final high-confidence viral sequences (both assemblers)
    
    ==========================================
    """
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Metagenome Viral Classification Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}

