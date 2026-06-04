process FRIP_SCORE {
    label 'process_low'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(peaks), path(flagstat) 
    
    output:
    tuple val(sample), path('*.frip.txt'), emit: frip
    
    script:
    """
    # Get total mapped reads from flagstat
    TOTAL=\$(grep "mapped (" ${flagstat} | head -1 | awk '{print \$1}')
    
    # Count reads in peaks - use -u flag to get unique reads that overlap
    IN_PEAKS=\$(bedtools intersect -a ${bam} -b ${peaks} -u -bed | wc -l)
    
    # Calculate FRiP using awk
    FRIP=\$(awk -v ip=\$IN_PEAKS -v total=\$TOTAL 'BEGIN {printf "%.4f", ip / total}')
    
    echo -e "${sample}\\t\$TOTAL\\t\$IN_PEAKS\\t\$FRIP" > ${sample}.frip.txt
    """
    
    stub:
    """
    touch ${sample}.frip.txt
    """
}