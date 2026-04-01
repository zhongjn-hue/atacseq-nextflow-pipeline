process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::macs2=2.1.2'
    container 'biocontainers/macs2:2.1.2--py39hf95cd2a_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.narrowPeak'),  emit: peaks
    tuple val(meta), path('*.summits.bed'), emit: summits
    tuple val(meta), path('*.xls'),         emit: xls
    path 'versions.yml',                    emit: versions

    script:
    def prefix = meta.id
    """
    macs2 callpeak \\
        -t ${bam} \\
        -f BAMPE \\
        -n ${prefix} \\
        --outdir . \\
        -g mm \\
        --nomodel \\
        --shift -100 \\
        --extsize 200 \\
        --nolambda \\
        --keep-dup all \\
        -B \\
        --SPMR \\
        --call-summits \\
        2>&1 | tee ${prefix}_macs2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version 2>&1 | sed 's/macs2 //')
    END_VERSIONS
    """
}