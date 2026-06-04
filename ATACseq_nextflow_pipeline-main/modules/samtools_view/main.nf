#!/usr/bin/env nextflow

process SAMTOOLS_VIEW {
    label 'process_low'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: "copy"
    
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path('*.filtered.bam'), emit: bam

    script:
    """
    samtools view -h -F 4 -q 30 ${bam} | \\
        grep -v chrM | \\
        samtools view -b -@ ${task.cpus} -o ${sample}.filtered.bam
    """

    stub:
    """
    touch ${sample}.filtered.bam
    """
}