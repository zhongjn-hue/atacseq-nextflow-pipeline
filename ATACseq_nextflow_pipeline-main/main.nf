include { SRA_DOWNLOAD } from './modules/sra_tools'
include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { BOWTIE2_BUILD } from './modules/bowtie2_build'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat'
include { SAMTOOLS_SORT } from './modules/samtools_sort'
include { SAMTOOLS_VIEW } from './modules/samtools_view'
include { SAMTOOLS_IDX } from './modules/samtools_idx'
include { MULTIQC } from './modules/multiqc'
include { TRIM } from './modules/trimmomatic'
include { TAGDIR } from './modules/homer_maketagdir'
include { FINDPEAKS } from './modules/homer_findpeaks'
include { POS2BED } from './modules/homer_pos2bed'
include { BEDTOOLS_MERGE } from './modules/bedtools_merge'
include { BEDTOOLS_INTERSECT } from './modules/bedtools_intersect'
include { FRIP_SCORE } from './modules/frip_score'
include { TSS_ENRICHMENT } from './modules/tss_enrich'
include { SUMMARIZE_TSS } from './modules/tss_summary'
include { BAMCOVERAGE } from './modules/deeptools_bamcoverage'
include { COMPUTEMATRIX } from './modules/deeptools_computematrix'


workflow {

    // Read samples CSV and create tuples (sample, srr)
    Channel.fromPath(params.sra)
        | splitCsv(header: true)
        | map { row -> tuple(row.sample, row.srr) }
        | set { download_ch }

    // Download FASTQ files
    downloaded_fastq = SRA_DOWNLOAD(download_ch)

    // Run FastQC on downloaded FASTQ
    FASTQC_RAW(downloaded_fastq)

    // Trimming off the adapters I found in FastQC
    trimmed_fastq = TRIM(downloaded_fastq, params.adapters)

    // Run FastQC on trimmed FASTQ to verify adapter removal
    FASTQC_TRIMMED(trimmed_fastq.trimmed_reads)

    // Create the Bowtie Index -> used same genome .fa as paper
    BOWTIE2_BUILD(params.genome)

    // Aligning the reads to the genome - use TRIMMED reads now
    BOWTIE2_ALIGN(trimmed_fastq.trimmed_reads, BOWTIE2_BUILD.out.index, BOWTIE2_BUILD.out.name) 

    // Getting the stats for the quality of the alignment
    SAMTOOLS_FLAGSTAT(BOWTIE2_ALIGN.out.bam) 

    // Make the channel to run to MultiQC
    multiqc_ch = FASTQC_RAW.out.zip
    .mix(FASTQC_TRIMMED.out.zip)
    .mix(trimmed_fastq.log)
    .mix(SAMTOOLS_FLAGSTAT.out.flagstat)
    .map { tuple -> tuple[1] }          // extract the file from the tuple
    .collect()
    .map { files -> files.flatten() }   // flatten nested lists into a single list

    // MultiQC to check the quality first
    MULTIQC(multiqc_ch)

    // Sort the BAM
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)

    // Keep only high quality reads, and remove mito reads
    SAMTOOLS_VIEW(SAMTOOLS_SORT.out.bam)

    // Index the BAM
    SAMTOOLS_IDX(SAMTOOLS_VIEW.out.bam)

    // Make Homer Tag Directory from filtered BAM
    TAGDIR(SAMTOOLS_VIEW.out.bam)  // Use filtered BAM!
    tagdir_ch = TAGDIR.out.tagdir
    //tagdir_ch.view()

    // Find peaks with ATAC-seq style
    FINDPEAKS(tagdir_ch)

    // Convert to BED format
    bed_ch = POS2BED(FINDPEAKS.out)
        .map { sample, bed -> bed }
        .collect()
    //bed_ch.view()

    // Create consensus peaks
    BEDTOOLS_MERGE(bed_ch)
    //BEDTOOLS_MERGE.out.bed.view { "Consensus peaks file: $it" }

    // Count reads in consensus peaks for each sample
    BEDTOOLS_INTERSECT(SAMTOOLS_VIEW.out.bam, BEDTOOLS_MERGE.out.bed)

    // QC Metrics: FRiP Score
    frip_input = SAMTOOLS_VIEW.out.bam
    .join(POS2BED.out)
    .join(SAMTOOLS_FLAGSTAT.out.flagstat)
    // frip_input.view()

    FRIP_SCORE(frip_input)

    //QC Metrics: TSS Enrichment
    gtf = file(params.gtf)
    TSS_ENRICHMENT(TAGDIR.out.tagdir,gtf)
    SUMMARIZE_TSS(TSS_ENRICHMENT.out.enrichment.map { it[1] }.collect())

    //Trying to turn the BAM --> BigWig
    BAMCOVERAGE(SAMTOOLS_IDX.out.index)

    // Compute TSS enrichment matrix for each sample individually
    tss_bed = file(params.ucsc_genes)
    COMPUTEMATRIX(BAMCOVERAGE.out.bigwig, tss_bed)

    // Your  COMPUTEMATRIX runs per sample
    matrices = COMPUTEMATRIX.out  // tuple(sample_id, matrix.gz)
    
}