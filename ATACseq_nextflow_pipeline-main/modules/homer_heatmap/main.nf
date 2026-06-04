#!/usr/bin/env nextflow

process HOMER_HEATMAP_MATRIX {
    label 'process_medium'
    publishDir "${params.outdir}/homer_heatmaps", mode: 'copy', pattern: "*.txt"
    container 'ghcr.io/bf528/homer_samtools:latest'
    stageInMode 'copy'
    
    input:
    tuple val(cell_type), path(gain_bed), path(loss_bed), path('tagdirs/*')
    
    output:
    tuple val(cell_type), path("*_matrix.txt"), emit: matrices
    
    script:
    """
    # Create matrices for each tag directory
    for tagdir in tagdirs/*; do
        sample=\$(basename \$tagdir | sed 's/_tags//')
        
        echo "Processing \$sample for ${cell_type}..."
        
        annotatePeaks.pl ${gain_bed} ${params.genome} \\
            -size 3000 -hist 25 \\
            -d \$tagdir \\
            > \${sample}_gain_matrix.txt
        
        annotatePeaks.pl ${loss_bed} ${params.genome} \\
            -size 3000 -hist 25 \\
            -d \$tagdir \\
            > \${sample}_loss_matrix.txt
    done
    
    echo "Matrix generation complete for ${cell_type}"
    """
}