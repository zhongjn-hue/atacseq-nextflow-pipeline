process BOWTIE2_INDEX {
    tag "$fasta"
    label 'process_high'

    conda 'bioconda::bowtie2=2.5.4'
    container 'biocontainers/bowtie2:2.5.4--he20e202_1'

    input:
    path fasta

    output:
    path 'bowtie2_index', emit: index
    path 'versions.yml',  emit: versions

    script:
    """
    mkdir bowtie2_index
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${fasta} \\
        bowtie2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*bowtie2-align-s version //')
    END_VERSIONS
    """
}