#!/usr/bin/env nextflow

process BEDTOOLS_INTERSECT {
    label 'process_single'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(bam)
    path(consensus_peaks)

    output:
    tuple val(sample), path('*.counts.txt'), emit: counts
    path('*.debug.log'), emit: log

    script:
    """
    echo "Consensus peaks info:" > ${sample}.debug.log
    head -5 ${consensus_peaks} >> ${sample}.debug.log
    wc -l ${consensus_peaks} >> ${sample}.debug.log
    
    echo "BAM info:" >> ${sample}.debug.log
    samtools view -H ${bam} | grep '@SQ' | head -5 >> ${sample}.debug.log
    samtools view ${bam} | head -5 >> ${sample}.debug.log
    
    echo "Running bedtools intersect..." >> ${sample}.debug.log
    bedtools intersect \\
        -a ${consensus_peaks} \\
        -b ${bam} \\
        -c > ${sample}.counts.txt
    
    echo "Output info:" >> ${sample}.debug.log
    wc -l ${sample}.counts.txt >> ${sample}.debug.log
    head -5 ${sample}.counts.txt >> ${sample}.debug.log
    """

    stub:
    """
    touch ${sample}.counts.txt
    touch ${sample}.debug.log
    """
}