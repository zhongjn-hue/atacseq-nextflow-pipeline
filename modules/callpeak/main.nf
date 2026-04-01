process BEDTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::bedtools=2.31.1'
    container 'biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    input:
    tuple val(meta), path(peaks)
    path blacklist   // mm10 ENCODE blacklist BED

    output:
    tuple val(meta), path('*_filtered.narrowPeak'), emit: filtered_peaks
    tuple val(meta), path('*_peak_stats.txt'),       emit: stats
    path 'versions.yml',                             emit: versions

    script:
    def prefix = meta.id
    """
    bedtools intersect \\
        -v \\
        -a ${peaks} \\
        -b ${blacklist} \\
        > ${prefix}_filtered.narrowPeak

    TOTAL=\$(wc -l < ${peaks})
    FILTERED=\$(wc -l < ${prefix}_filtered.narrowPeak)
    REMOVED=\$((TOTAL - FILTERED))
    echo "Total peaks before filtering: \$TOTAL"    > ${prefix}_peak_stats.txt
    echo "Peaks after blacklist removal: \$FILTERED" >> ${prefix}_peak_stats.txt
    echo "Blacklisted peaks removed: \$REMOVED"     >> ${prefix}_peak_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}