#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'
    container 'ghcr.io/bf528/multiqc:latest'
    publishDir params.outdir, mode: "copy"

    input:
    path('*')

    output:
    path('*.html'), emit: html

    script:
    """
    multiqc . -f --outdir . --filename multiqc_report.html
    """

    stub:
    """
    touch multiqc_report.html
    """
}

