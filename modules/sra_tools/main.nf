#!/usr/bin/env nextflow

process SRA_DOWNLOAD {

    label 'process_veryhigh'
    container 'ghcr.io/bf528/sratools:latest'

    input:
    tuple val(sample), val(srr)

    output:
    tuple val(sample), path("${srr}.fastq.gz")

    script:
    """
    mkdir -p ${params.data_dir}

    # Download FASTQ
    fasterq-dump ${srr} --threads ${task.cpus} -O ${params.data_dir}

    # Compress the downloaded FASTQ
    gzip ${params.data_dir}/${srr}.fastq

    # Symlink to work directory for Nextflow
    ln -s ${params.data_dir}/${srr}.fastq.gz .
    """

    stub:
    """
    touch ${srr}.fastq.gz
    """
}

