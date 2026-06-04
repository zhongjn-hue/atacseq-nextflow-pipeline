#!/usr/bin/env nextflow

process SUMMARIZE_TSS {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/tss_enrichment", mode: 'copy'
    
    input:
    path(enrichment_files)
    
    output:
    path("tss_enrichment_summary.tsv")
    
    script:
    """
    echo -e "Sample\\tTSS_Enrichment\\tStatus" > tss_enrichment_summary.tsv
    
    for file in *_tss_enrichment.txt; do
        sample=\$(grep "^Sample" \$file | cut -f2)
        enrich=\$(grep "^TSS_Enrichment" \$file | cut -f2)
        status=\$(grep "^Status" \$file | cut -f2)
        echo -e "\$sample\\t\$enrich\\t\$status" >> tss_enrichment_summary.tsv
    done
    """
}