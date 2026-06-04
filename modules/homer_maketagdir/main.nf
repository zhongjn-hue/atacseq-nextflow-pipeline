#!/usr/bin/env nextflow

process TAGDIR {
    label 'process_single'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode: "copy"

    input: 
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_tags"), emit: tagdir

    script:
    """
    makeTagDirectory ${sample_id}_tags ${bam_file}
    """

    stub:
    """
    mkdir -p ${sample_id}_tags
    echo "stub tag info for ${sample_id}" > ${sample_id}_tags/tagInfo.txt
    """
}
