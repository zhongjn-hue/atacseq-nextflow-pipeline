#!/usr/bin/env nextflow

process POS2BED {
    label 'process_single'
    container 'ghcr.io/bf528/homer:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(peaks)

    output:
    tuple val(sample), path('*.bed'), emit: bed

    script:
    """
    pos2bed.pl ${peaks} > ${sample}.bed
    """

    stub:
    """
    touch ${sample}.bed
    """
}