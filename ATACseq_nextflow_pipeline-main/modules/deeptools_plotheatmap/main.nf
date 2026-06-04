#!/usr/bin/env nextflow

process HEATMAP {
    label 'process_single'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(group_id), path('gain_matrix.gz'), path('loss_matrix.gz')

    output:
    path("${group_id}_differential_heatmap.png")

    script: 
    """
    # Relabel the groups so they have distinct names
    computeMatrixOperations relabel \
        -m gain_matrix.gz \
        -o gain_relabeled.gz \
        --groupLabels "Gained"
    
    computeMatrixOperations relabel \
        -m loss_matrix.gz \
        -o loss_relabeled.gz \
        --groupLabels "Lost"
    
    # Now combine them - they'll maintain separate groups
    computeMatrixOperations rbind \
        -m gain_relabeled.gz loss_relabeled.gz \
        -o ${group_id}_differential_matrix.gz
    
    # Plot combined heatmap with two distinct sections
    plotHeatmap \
        --matrixFile ${group_id}_differential_matrix.gz \
        --outFileName ${group_id}_differential_heatmap.png \
        --colorMap Reds \
        --missingDataColor white \
        --heatmapHeight 15 \
        --whatToShow 'heatmap and colorbar' \
        --sortRegions descend \
        --sortUsing mean \
        --refPointLabel "TSS" \
        --plotTitle "${group_id} Differential Chromatin Accessibility" \
        --xAxisLabel "Distance from TSS (bp)"
    """

    stub:
    """
    touch ${group_id}_differential_heatmap.png
    """
}