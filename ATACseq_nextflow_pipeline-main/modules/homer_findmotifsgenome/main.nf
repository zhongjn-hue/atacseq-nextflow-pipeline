#!/usr/bin/env nextflow

process FIND_MOTIFS_GENOME {
    label 'process_high'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode: "copy"

    input: 
    tuple val(cell_type), path(significant_peaks_bed)
    path(genome)

    output: 
    tuple val(cell_type), path("${cell_type}_motifs")

    script:
    """
    mkdir ${cell_type}_motifs

    findMotifsGenome.pl ${significant_peaks_bed} ${genome} ${cell_type}_motifs \
    -size 200 \
    -mask \
    -p ${task.cpus} \
    -len 8,10,12

    """

    stub:
    """
    mkdir ${cell_type}_motifs
    """
}