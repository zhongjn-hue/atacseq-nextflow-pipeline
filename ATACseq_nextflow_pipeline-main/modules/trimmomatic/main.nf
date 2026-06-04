process TRIM {
    label 'process_medium'
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample_id), path(read)
    path(adapters)

    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.log"), emit: log

    script:
    """
    trimmomatic SE -threads ${task.cpus} \\
        ${read} \\
        ${sample_id}_trimmed.fastq.gz \\
        ILLUMINACLIP:${adapters}:2:20:7 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \\
        2> ${sample_id}_trimmomatic.log
    """

    stub:
    """
    touch ${sample_id}_trimmed.fastq.gz
    touch ${sample_id}_trimmomatic.log
    """
}