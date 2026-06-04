//Run the R Script for DiffBind BEFORE this analysis

include { ANNOTATE } from './modules/homer_annotatepeaks'
include { FIND_MOTIFS_GENOME } from './modules/homer_findmotifsgenome'
include { COMPUTEMATRIX_DIFFERENTIAL } from './modules/deeptools_computematrixdiff'
include { HEATMAP } from './modules/deeptools_plotheatmap'

workflow {

  ch_significant_peaks = Channel.of(
      ['cDC1', file("${params.diffbind_dir}/cdc1_significant_peaks.bed")],
      ['cDC2', file("${params.diffbind_dir}/cdc2_significant_peaks.bed")]
  )

  // Channel for gain/loss BED files
  ch_differential_beds = Channel.of(
      ['cDC1', file("${params.diffbind_dir}/cdc1_gain_peaks.bed"), file("${params.diffbind_dir}/cdc1_loss_peaks.bed")],
      ['cDC2', file("${params.diffbind_dir}/cdc2_gain_peaks.bed"), file("${params.diffbind_dir}/cdc2_loss_peaks.bed")]
  )

  // Channel for bigwig files (from main.nf output)
  ch_bigwigs = Channel.fromPath("${params.outdir}/*.bw")
      .map { bw ->
          def sample_id = bw.name.replaceAll('.bw', '')
          def cell_type = sample_id.split('_')[0]
          tuple(cell_type, bw)
      }
      .groupTuple()

  // Join bigwigs with differential beds by group_id (cDC1/cDC2)
  ch_matrix_input = ch_bigwigs.join(ch_differential_beds)
  //ch_matrix_input.view()

  // Compute matrix for differential regions (gain + loss)
  COMPUTEMATRIX_DIFFERENTIAL(ch_matrix_input)
  //COMPUTEMATRIX_DIFFERENTIAL.out.view()

  // Generate differential accessibility heatmaps
  HEATMAP(COMPUTEMATRIX_DIFFERENTIAL.out)

  // Use .collect() to make it a value channel that can be reused
  genome_fasta = Channel.fromPath(params.genome).collect()
  gtf_file = Channel.fromPath(params.gtf).collect()
  
  // Find enriched motifs in significant peaks
  FIND_MOTIFS_GENOME(ch_significant_peaks, genome_fasta)

  // Annotate peaks with genomic features and gene associations
  ANNOTATE(ch_significant_peaks, genome_fasta, gtf_file)
  
}