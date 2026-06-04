#!/usr/bin/env nextflow

process COMPUTEMATRIX_DIFFERENTIAL {
    label 'process_veryhigh'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(group_id), path(bigwigs), path(gain_bed), path(loss_bed)

    output:
    tuple val(group_id), path("${group_id}_gain_matrix.gz"), path("${group_id}_loss_matrix.gz")

    script:
    def window = 2000
    """
    # Create matrix for GAIN regions
    computeMatrix reference-point \
        -S ${bigwigs} \
        -R ${gain_bed} \
        --referencePoint TSS \
        --beforeRegionStartLength ${window} \
        --afterRegionStartLength ${window} \
        --binSize 10 \
        --missingDataAsZero \
        --numberOfProcessors ${task.cpus} \
        -o ${group_id}_gain_matrix.gz
    
    # Create matrix for LOSS regions  
    computeMatrix reference-point \
        -S ${bigwigs} \
        -R ${loss_bed} \
        --referencePoint TSS \
        --beforeRegionStartLength ${window} \
        --afterRegionStartLength ${window} \
        --binSize 10 \
        --missingDataAsZero \
        --numberOfProcessors ${task.cpus} \
        -o ${group_id}_loss_matrix.gz
    """

    stub:
    """
    touch ${group_id}_gain_matrix.gz
    touch ${group_id}_loss_matrix.gz
    """
}