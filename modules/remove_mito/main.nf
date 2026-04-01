#!/usr/bin/env nextflow
process REMOVE_MITO {
    label 'process_medium'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/04_filtered_bam", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), val(celltype), path(bam)

    output:
    tuple val(sample_id), val(condition), val(celltype), path("${sample_id}_filtered.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_filter_stats.txt"), emit: stats

    script:
    """
    samtools view -h ${bam} | grep -v chrM | samtools view -b > ${sample_id}_filtered.bam
    samtools index ${sample_id}_filtered.bam
    
    echo "Original reads: \$(samtools view -c ${bam})" > ${sample_id}_filter_stats.txt
    echo "Filtered reads: \$(samtools view -c ${sample_id}_filtered.bam)" >> ${sample_id}_filter_stats.txt
    """
}