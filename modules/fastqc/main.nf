#!/usr/bin/env nextflow

process FASTQC {

    label 'process_low'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir "${params.outdir}/fastqc", mode: "copy"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*.zip'), emit: zip
    tuple val(sample), path('*.html'), emit: html

    script:
    """
    export _JAVA_OPTIONS="-Djava.awt.headless=true"
    fastqc $reads -t $task.cpus
    """

    stub:
    """
    touch stub_${sample}_fastqc.zip
    touch stub_${sample}_fastqc.html
    """

}