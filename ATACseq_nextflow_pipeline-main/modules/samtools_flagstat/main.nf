#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    label 'process_single'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample), path(bam_file)

    output:
    tuple val(sample), path("${sample}.flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat ${bam_file} > ${sample}.flagstat.txt
    """

    stub:
    """
    touch ${sample}.flagstat.txt
    """
}
