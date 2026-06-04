#!/usr/bin/env nextflow

process COMPUTEMATRIX {
    label 'process_veryhigh'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(sample_id), path(bigwig)
    path(bed)

    output:
    tuple val(sample_id), path("${sample_id}_tss_matrix.gz")

    script:
    def window = 2000
    """
    computeMatrix reference-point \
        -S ${bigwig} \
        -R ${bed} \
        --referencePoint TSS \
        --beforeRegionStartLength ${window} \
        --afterRegionStartLength ${window} \
        --binSize 10 \
        --missingDataAsZero \
        --numberOfProcessors ${task.cpus} \
        -o ${sample_id}_tss_matrix.gz
    """

    stub:
    """
    touch ${sample_id}_tss_matrix.gz
    """
}