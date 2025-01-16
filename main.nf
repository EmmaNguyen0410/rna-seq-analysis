#!/usr/bin/env nextflow

params.outdir = "results"
params.fastq = "$baseDir/data/reads/*/*_{R1,R2}_*.fastq.gz"
params.ref = "$baseDir/data/GCA_027886375.1_ASM2788637v1_genomic.fna"
params.gtf = "$baseDir/data/GCA_027886375.1_ASM2788637v1_genomic.gtf"
params.bowtie2_index = "$baseDir/data/bowtie2_index"
params.script = "$baseDir/bin/deseq2.R"


process FASTQC {
    executor 'local'

    // resources_used.mem = 1205420kb
    // resources_used.ncpus = 2
    memory '2 GB'
    cpus 2
    tag "FASTQC on $sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "$params.outdir/fastqc", mode: 'move'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs", emit: logs

    script:
    """
    echo "NXF_TASK_WORKDIR: \$#NXF_TASK_WORKDIR"
    echo "TMPDIR: \$#TMPDIR"
    mkdir fastqc_${sample_id}_logs
    fastqc -t 2 -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process TRIMMING {
    executor 'pbspro'

    // resources_used.mem = 356704kb
    // resources_used.ncpus = 8
    memory '1 GB'
    cpus 8
    
    tag "Trimming $sample_id"
    conda 'bioconda::cutadapt=4.9'
    // publishDir "$params.outdir/trimming", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_{1,2}.fastq.gz"), emit: trimmed

    script:
    """
    cutadapt --compression-level=2 -q 20 -m 20 -j 8 -o trimmed_${sample_id}_1.fastq.gz -p trimmed_${sample_id}_2.fastq.gz $reads
    """
}

process MAPPING {
    executor 'pbspro'

    // resources_used.mem = 278212kb
    // resources_used.ncpus = 9
    memory '1 GB'
    cpus 9
    
    tag "Mapping $sample_id"
    conda 'bioconda::bowtie2=2.5.4 bioconda::samtools=1.21'
    // publishDir "$params.outdir/mapped", mode: 'copy'

    input:
    path reference_genome
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("mapped_${sample_id}.bam"), emit: bam

    script:
    """
    bowtie2 -p 8 --very-sensitive -x $params.bowtie2_index/gonorrhea -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} -S mapped_${sample_id}.sam
    samtools view -S -b mapped_${sample_id}.sam > mapped_${sample_id}.bam

    rm ${trimmed_reads[0]} \$(readlink -f ${trimmed_reads[0]})
    rm ${trimmed_reads[1]} \$(readlink -f ${trimmed_reads[1]})
    rm mapped_${sample_id}.sam
    """
}

process QUANTIFICATION {
    executor 'pbspro'

    // resources_used.mem = 959796kb
    // resources_used.ncpus = 1
    memory '1 GB' 
    
    tag "Quantifying $sample_id"
    conda 'bioconda::htseq=2.0.5'
    publishDir "$params.outdir/quantification", mode: 'copy'

    input:
    path gtf
    tuple val(sample_id), path(bam)

    output:
    path "counts_${sample_id}.txt", emit: counts

    script:
    """
    samtools sort -n -o mapped_${sample_id}_sorted.bam ${bam} 
    htseq-count -f bam -r name -m union -s reverse -t gene -i gene_id mapped_${sample_id}_sorted.bam ${gtf} > counts_${sample_id}.txt
    
    rm ${bam} \$(readlink -f ${bam})
    rm mapped_${sample_id}_sorted.bam
    """
}

process DIFFERENTIAL_EXPRESSION {
    executor 'local'

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

    FASTQC(read_pairs_ch)
    TRIMMING(read_pairs_ch)
    MAPPING(params.ref, TRIMMING.out)    
    QUANTIFICATION(params.gtf, MAPPING.out)
    
    QUANTIFICATION.out
    .collect()
    .set{ htseq_outputs }

    DIFFERENTIAL_EXPRESSION(htseq_outputs)
}