#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    label 'process_veryhigh'
    container 'ghcr.io/bf528/bowtie2:latest'

    input:
    tuple val(sample), path(read)
    path bt2
    val name

    output:
    tuple val(sample), path("*.bam"), emit: bam

    shell:
    """
    # Copy index files into current working directory so bowtie2 can find them
    cp -r ${bt2}/* .

    # Run alignment (single-end)
    bowtie2 -p 8 --very-sensitive -x ${name} -U ${read} | samtools view -bS - > ${sample}.bam
    """

    stub:
    """
    touch ${sample}.bam
    """
}