#!/usr/bin/env nextflow

process SAMTOOLS_IDX {
    label 'process_single'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(sorted_bam)

    output:
    tuple val(sample), path(sorted_bam), path("${sorted_bam}.bai"), emit: index

    script:
    """
    samtools index ${sorted_bam}
    """
    
    stub:
    """
    touch ${sorted_bam}.bai
    """
}

