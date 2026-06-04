#!/usr/bin/env nextflow

process FASTQC {
    label 'process_low'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir params.outdir, mode: "copy"
    
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*.html'), emit: html
    tuple val(sample), path('*.zip'), emit: zip

    script:
    """
    fastqc -t $task.cpus $reads
    """

    stub:
    """
    touch ${sample}_fastqc.html
    touch ${sample}_fastqc.zip
    """

}

