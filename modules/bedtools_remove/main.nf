#!/usr/bin/env nextflow

process BEDTOOLS_REMOVE {
    label 'process_low'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path(repr_peaks_bed)
    path(blacklist)

    output:
    path('repr_peaks_filtered.bed'), emit: bed

    script:
    """
    bedtools intersect -a ${repr_peaks_bed} -b ${blacklist} -v > repr_peaks_filtered.bed
    """   

    stub:
    """
    touch repr_peaks_filtered.bed
    """
}