#!/usr/bin/env nextflow

params.outdir = "results"
params.fastq = "$baseDir/data/test/*_{1,2}.fq"
// params.fastq = "$baseDir/data/test/*.fq"
params.multiqc = "$baseDir/multiqc"
params.ref = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.gff = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.script = "$baseDir/bin/deseq2.R"

// process FASTQC {
//     memory '32 GB'
//     cpus 4
//     tag "FASTQC on $sample_id"
//     conda 'bioconda::fastqc=0.12.1'
//     publishDir "$params.outdir/fastqc", mode: 'copy'

//     input:
//     tuple val(sample_id), path(reads)

//     output:
//     path "fastqc_${sample_id}_logs", emit: logs

//     script:
//     """
//     mkdir fastqc_${sample_id}_logs
//     fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
//     """
// }

process FASTQC {
    memory '4 GB'
    cpus 2
    tag "FASTQC on"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "$params.outdir/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs", emit: logs

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t 2 -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

// process MULTIQC {
//     tag "MultiQC $sample_id"
//     conda 'bioconda::multiqc=1.25'
//     publishDir "$params.outdir/multiqc", mode:'copy'

//     input: 
//     path '*'
//     path config

//     output:
//     path 'multiqc_report.html', emit: report

//     script:
//     """
//     echo "custom_logo: ./multiqc/logo.png" >> ./multiqc/multiqc_config.yaml
//     multiqc -o multiqc_report.html .
//     """
// }

process TRIMMING {
    memory '32 GB'
    cpus 8
    tag "Trimming $sample_id"
    conda 'bioconda::cutadapt=4.9'
    publishDir "$params.outdir/trimming", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_{1,2}.fastq"), emit: trimmed

    script:
    """
    cutadapt -q 20 -m 20 -j 8 -o trimmed_${sample_id}_1.fastq -p trimmed_${sample_id}_2.fastq $reads
    """
}

process MAPPING {
    memory '32 GB'
    cpus 8
    tag "Mapping $sample_id"
    conda 'bioconda::bowtie2=2.5.4 bioconda::samtools=1.21'
    publishDir "$params.outdir/mapped", mode: 'copy'

    input:
    path reference_genome
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("mapped_${sample_id}.bam"), emit: bam

    script:
    """
    mkdir bowtie2_index
    cd bowtie2_index
    bowtie2-build ../${reference_genome} ggal 
    cd .. 
    bowtie2 -p 8 --very-sensitive -x bowtie2_index/ggal -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} -S mapped_${sample_id}.sam
    samtools view -S -b mapped_${sample_id}.sam > mapped_${sample_id}.bam
    """
}

process QUANTIFICATION {
    memory '32 GB'
    cpus 8    
    tag "Quantifying $sample_id"
    conda 'bioconda::htseq=2.0.5'
    publishDir "$params.outdir/quantification", mode: 'copy'

    input:
    path gff
    tuple val(sample_id), path(bam)

    output:
    path "counts_${sample_id}.txt", emit: counts

    script:
    """
    samtools sort -n -o mapped_${sample_id}_sorted.bam ${bam} 
    htseq-count -n 3 -f bam -r name -m union -s reverse -t exon -i gene_id mapped_${sample_id}_sorted.bam ${gff} > counts_${sample_id}.txt
    """
}

process DIFFERENTIAL_EXPRESSION {
    memory '32 GB'
    cpus 3
    tag "Differential Expression"
    conda 'bioconda::deseq2=1.42.0'
    publishDir "$params.outdir/deseq2", mode: 'copy'

    input:
    path htseq_outputs 

    output:
    path "deseq2.pdf", emit: deseq2_results

    script:
    """
    mkdir htseq_outputs
    mv counts_*.txt ./htseq_outputs
    Rscript ${params.script} './htseq_outputs' 
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.fastq, checkIfExists: true)
    //read_pairs_ch = Channel.fromPath(params.fastq)
    FASTQC(read_pairs_ch)
    //MULTIQC(FASTQC.out, params.multiqc)
    // TRIMMING(read_pairs_ch)
    // MAPPING(params.ref, TRIMMING.out)
    // QUANTIFICATION(params.gff, MAPPING.out)

    // QUANTIFICATION.out
    // .collect()
    // .set{ htseq_outputs }

    // htseq_outputs = ['/shared/homes/162557/home/uts/speedx/RNA-seq-2/data/quantification/counts_ggal_gut_treated.txt', '/shared/homes/162557/home/uts/speedx/RNA-seq-2/data/quantification/counts_ggal_liver_untreated.txt', '/shared/homes/162557/home/uts/speedx/RNA-seq-2/data/quantification/counts_ggal_lung_treated.txt', '/shared/homes/162557/home/uts/speedx/RNA-seq-2/data/quantification/counts_ggal_spleen_untreated.txt']

    // DIFFERENTIAL_EXPRESSION(htseq_outputs)
}