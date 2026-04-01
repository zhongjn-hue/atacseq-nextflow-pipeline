process BAM_TO_BIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::deeptools=3.5.4'
    container 'biocontainers/deeptools:3.5.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.bw'), emit: bigwig
    path 'versions.yml',           emit: versions

    script:
    def prefix = meta.id
    // mm10 effective genome size = 2,652,783,500
    """
    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${prefix}.bw \\
        --outFileFormat bigwig \\
        --binSize 10 \\
        --normalizeUsing RPGC \\
        --effectiveGenomeSize 2652783500 \\
        --ignoreForNormalization chrX chrM \\
        --extendReads \\
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed 's/bamCoverage //')
    END_VERSIONS
    """
}