# RNA seq pipeline with Nextflow, Apptainer, and PBS 

## Tools and their pararmeters configuration 

This is a list of all tools used in this project and explanation for their usage: 

1. Fastqc

```
fastqc -t 2 -o fastqc_${sample_id}_logs -f fastq -q ${reads}
```
__-t__: works on file basis i.e one thread per file. 

2. [Cutadapt](https://cutadapt.readthedocs.io/en/v4.8/guide.html#basic-usage)

```
cutadapt -q 20 -m 20 -j 8 -o trimmed_${sample_id}_1.fastq -p trimmed_${sample_id}_2.fastq $reads
```

Whether an input file needs to be decompressed or an output file needs to be compressed is detected automatically by inspecting the file name. Because all output files are short-lived intermediate files, so they are not compressed to speed up the process (output file not ending with .gz). 

__-q__ (or --quality-cutoff): trim low-quality ends from reads. If you specify a single cutoff value, the 3’ end of each read is trimmed. For Illumina reads, this is sufficient as their quality is high at the beginning, but degrades towards the 3’ end.

__-m__: Discard processed reads that are shorter than LENGTH.

__-j 8__: Increasing the number of cores with -j will increase the number of reads per minute at near-linear rate.

__-p__: By default, all processed reads, no matter whether they were trimmed or not, are written to the output file specified by the -o option (or to standard output if -o was not provided). For paired-end reads, the second read in a pair is always written to the file specified by the -p option.

3. [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

4. [Samtools](https://www.htslib.org/doc/samtools-sort.html)

```
samtools sort -n -o mapped_${sample_id}_sorted.bam ${bam} 
```

__-n__: Sort by read names (i.e., the QNAME field) using an alpha-numeric ordering, rather than by chromosomal coordinates. 

5. [Htseq](https://htseq.readthedocs.io/en/release_0.11.1/count.html)

```
htseq-count -n 2 -f bam -r name -m union -s reverse -t exon -i gene_id mapped_${sample_id}_sorted.bam ${gff} > counts_${sample_id}.txt
```

__-r name__: For paired-end data, the alignment have to be sorted either by read name or by alignment position. 

__-s reverse__: The second read has to be on the same strand and the first read on the opposite strand. 

6. [Deseq2](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
