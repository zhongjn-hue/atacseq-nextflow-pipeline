#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
    label 'process_single'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(bam_file)
    
    output:
    tuple val(sample), path("${sample}.sorted.bam"), emit: bam

    script:
    """
    samtools sort -o ${sample}.sorted.bam ${bam_file}
    """

    stub:
    """
    touch ${sample}.sorted.bam
    """
}
