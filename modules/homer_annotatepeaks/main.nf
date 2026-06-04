#!/usr/bin/env nextflow

process ANNOTATE {
    label 'process_high'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/annotation", mode: "copy"

    // The bed is the diffbind results
    input:
    tuple val(cell_type), path(bed)
    path(genome)
    path(gtf)

    output:
    tuple val(cell_type), path("${cell_type}_annotated_peaks.txt")

    script:
    """
    annotatePeaks.pl ${bed} ${genome} -gtf ${gtf} > ${cell_type}_annotated_peaks.txt
    """

    stub:
    """
    touch ${cell_type}_annotated_peaks.txt
    """
}