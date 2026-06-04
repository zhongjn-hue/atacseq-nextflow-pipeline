#!/usr/bin/env nextflow

process FINDPEAKS {
    label 'process_high'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode: "copy"

    input: 
    tuple val(sample_id), path(tagdir)

    output: 
    tuple val(sample_id), path("${sample_id}_peaks.txt"), emit: peaks

    script: 
    """
    findPeaks ${tagdir} -style histone -o ${sample_id}_peaks.txt
    """

    stub:
    """
    touch ${sample_id}_peaks.txt
    """
}