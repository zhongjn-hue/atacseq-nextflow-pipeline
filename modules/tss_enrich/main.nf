#!/usr/bin/env nextflow

process TSS_ENRICHMENT {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/tss_enrichment", mode: 'copy'

    input:
    tuple val(sample), path(tagdir)
    path(gtf)

    output:
    tuple val(sample), path("${sample}_tss_enrichment.txt"), emit: enrichment
    path("${sample}_tss_distribution.txt"), emit: distribution

    script:
    """
    # Extract TSS positions from GTF (first base of each gene)
    awk '\$3 == "gene" {
        if (\$7 == "+") {
            print \$1 "\\t" \$4 "\\t" \$4+1 "\\t" \$10 "\\t0\\t" \$7
        } else {
            print \$1 "\\t" \$5-1 "\\t" \$5 "\\t" \$10 "\\t0\\t" \$7
        }
    }' ${gtf} | sed 's/"//g' | sed 's/;//g' > tss_positions.bed
    
    echo "Created \$(wc -l tss_positions.bed | cut -f1 -d' ') TSS positions" >&2
    
    # Calculate TSS enrichment histogram using the BED file
    annotatePeaks.pl tss_positions.bed none \\
        -size 4000 \\
        -hist 10 \\
        -d ${tagdir} \\
        > ${sample}_tss_distribution.txt
    
    # Check if output was created
    if [ ! -s ${sample}_tss_distribution.txt ]; then
        echo "ERROR: Distribution file is empty!" >&2
        echo "Sample\\t${sample}" > ${sample}_tss_enrichment.txt
        echo "TSS_Enrichment\\tNA" >> ${sample}_tss_enrichment.txt
        echo "Status\\tERROR" >> ${sample}_tss_enrichment.txt
        exit 0
    fi
    
    # Calculate enrichment
    cat > calc_tss.awk <<'AWK_SCRIPT'
BEGIN {
    tss_sum=0; tss_count=0
    bg_sum=0; bg_count=0
}
/^#/ { next }
/^Distance/ { next }
{
    dist = \$1 + 0
    cov = \$2 + 0
    
    if (dist >= -100 && dist <= 100) {
        tss_sum += cov
        tss_count++
    }
    if ((dist >= -2000 && dist <= -1000) || (dist >= 1000 && dist <= 2000)) {
        bg_sum += cov
        bg_count++
    }
}
END {
    if (tss_count > 0 && bg_count > 0) {
        tss_mean = tss_sum / tss_count
        bg_mean = bg_sum / bg_count
        enrichment = (bg_mean > 0) ? tss_mean / bg_mean : 0
        
        if (enrichment >= 7) status = "EXCELLENT"
        else if (enrichment >= 5) status = "GOOD"
        else if (enrichment >= 3) status = "ACCEPTABLE"
        else status = "POOR"
        
        print "Sample\\t${sample}"
        print "TSS_Enrichment\\t" enrichment
        print "Status\\t" status
        print "TSS_Coverage\\t" tss_mean
        print "Background_Coverage\\t" bg_mean
    } else {
        print "Sample\\t${sample}"
        print "TSS_Enrichment\\tNA"
        print "Status\\tERROR"
    }
}
AWK_SCRIPT

    awk -f calc_tss.awk ${sample}_tss_distribution.txt > ${sample}_tss_enrichment.txt
    
    cat ${sample}_tss_enrichment.txt
    """
}