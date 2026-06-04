#!/usr/bin/env nextflow

process BEDTOOLS_MERGE {
    label 'process_medium'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path(bed_files)

    output:
    path('consensus_peaks.bed'), emit: bed
    path('debug.log'), emit: log

    script:
    """
    # Extract first 3 columns and combine all files
    cat ${bed_files} | \\
        grep -v '^#' | \\
        awk '{print \$1"\t"\$2"\t"\$3}' > combined.bed
    
    echo "Combined bed lines:" > debug.log
    wc -l combined.bed >> debug.log
    
    echo "First 10 lines of combined:" >> debug.log
    head combined.bed >> debug.log
    
    # Sort and merge
    sort -k1,1 -k2,2n combined.bed | \\
        bedtools merge -i - > consensus_peaks.bed
    
    echo "Consensus peaks:" >> debug.log
    wc -l consensus_peaks.bed >> debug.log
    head consensus_peaks.bed >> debug.log
    """

    stub:
    """
    touch consensus_peaks.bed
    touch debug.log
    """
}