#!/usr/bin/env nextflow

process PLOTPROFILE {
    label 'process_low'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: "copy"

    input: 
    tuple val(sample_id), path(matrix)
    
    output: 
    path("${sample_id}_signal_coverage.png")

    script: 
    """
    plotProfile \
        -m ${matrix} \
        -out ${sample_id}_signal_coverage.png \
        --plotTitle "${sample_id} Signal Coverage" \
        --refPointLabel "TSS" \
        --regionsLabel "${sample_id}" \
        --perGroup
    """

    stub:
    """
    touch ${sample_id}_signal_coverage.png
    """
}