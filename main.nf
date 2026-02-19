include { FASTQC } from './modules/fastqc'
include { TRIMMOMATIC } from './modules/trimming'
include { BOWTIE2_INDEX } from './modules/bowtie2_index'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align'
include { REMOVE_MITO } from './modules/remove_mito'
include { SAMTOOLS_SORT } from './modules/samtools_sort'
include { SAMTOOLS_INDEX } from './modules/samtools_index'
include { MACS2_CALLPEAK } from './modules/marc2'
include { BEDTOOLS } from './modules/callpeak'
include { BAM_TO_BIGWIG } from './modules/bigwig'
include { MULTIQC } from './modules/multiqc'

// Parameters
params.samplesheet = 'samplesheet.csv'
params.genome_fasta = ''  // Will use value from nextflow.config
params.gtf = ''  // Will use value from nextflow.config
params.outdir = './results'
params.help = false

workflow {

    /* Build bowtie2 index */
    BOWTIE2_INDEX(file(params.genome_fasta))

    /* Read samplesheet */
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq)
            )
        }
        .set { samples_ch }

    /* Step 1: QC raw reads */
    FASTQC(samples_ch)

    /* Step 2: trimming */
    TRIMMOMATIC(samples_ch)

    /* Step 3: alignment */
    BOWTIE2_ALIGN(
        TRIMMOMATIC.out.trimmed,
        BOWTIE2_INDEX.out.index
    )

    /* Step 4: remove MT */
    REMOVE_MITO(BOWTIE2_ALIGN.out.bam)
    
    SAMTOOLS_SORT(REMOVE_MITO.out.bam)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    sorted_with_index = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    MACS2_CALLPEAK(sorted_with_index)

    /* Step 6: peak calling */
    BEDTOOLS(sorted_with_index)

    /* Step 7: bigwig */
    BAM_TO_BIGWIG(sorted_with_index)

    /* Final multiqc */
    qc_channel = FASTQC.out.html
                    .map { it[1] }
                    .mix(FASTQC.out.zip.map { it[1] })
                    .mix(BOWTIE2_ALIGN.out.stats.map { it[1] })
                    .mix(REMOVE_MITO.out.stats.map { it[1] })
                    .collect()
    
    MULTIQC(qc_channel)
}