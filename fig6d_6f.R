#!/usr/bin/env Rscript
library(Gviz)
library(rtracklayer)
library(GenomicRanges)

# Allow non-standard chromosome names
options(ucscChromosomeNames=FALSE)

# =============================================================================
# Configuration
# =============================================================================
# Paths
results_dir <- "results"
output_dir <- "results"
gtf_file <- "data/annot/gencode.vM10.primary_assembly.annotation.gtf"

# Color scheme
color_wt <- "black"
color_ko <- "#8B2323"  # Dark red/maroon
color_highlight_gray <- "#D3D3D3"
color_highlight_pink <- "#FFB6C1"

# =============================================================================
# Function to create browser tracks
# =============================================================================
create_browser_plot <- function(chr, start, end, gene_name, output_file, include_wt = FALSE) {
  
  cat(paste0("Creating plot for ", gene_name, " (", chr, ":", start, "-", end, ")\n"))
  
  # 1. Ideogram track (chromosome)
  itrack <- IdeogramTrack(genome = "mm10", chromosome = chr)
  
  # 2. Genome axis
  gtrack <- GenomeAxisTrack(
    col = "black",
    fontcolor = "black",
    fontsize = 10
  )
  
  # 3. Gene annotation track
  cat("Loading gene annotations...\n")
  grtrack <- GeneRegionTrack(
    gtf_file,
    genome = "mm10",
    chromosome = chr,
    start = start,
    end = end,
    name = "Genes",
    transcriptAnnotation = "symbol",
    background.title = "white",
    col = NULL,
    fill = "navy",
    fontcolor.group = "black",
    fontsize.group = 10
  )
  
  # 4. ATAC-seq data tracks
  cat("Loading bigWig files...\n")
  
  # Set y-axis limit based on gene
  y_max <- ifelse(gene_name == "Maged1", 150, 100)
  
  track_list <- list(itrack, gtrack, grtrack)
  
  # All samples, grouped by replicate (same order for both genes now)
  # cDC1 WT 1
  atac_cdc1_wt_1 <- DataTrack(
    range = file.path(results_dir, "cDC1_WT_1.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC1_WT_1",
    fill = color_wt,
    col = color_wt,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC1 KO 1
  atac_cdc1_ko_1 <- DataTrack(
    range = file.path(results_dir, "cDC1_KO_1.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC1_KO_1",
    fill = color_ko,
    col = color_ko,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC1 WT 2
  atac_cdc1_wt_2 <- DataTrack(
    range = file.path(results_dir, "cDC1_WT_2.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC1_WT_2",
    fill = color_wt,
    col = color_wt,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC1 KO 2
  atac_cdc1_ko_2 <- DataTrack(
    range = file.path(results_dir, "cDC1_KO_2.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC1_KO_2",
    fill = color_ko,
    col = color_ko,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC2 WT 1
  atac_cdc2_wt_1 <- DataTrack(
    range = file.path(results_dir, "cDC2_WT_1.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC2_WT_1",
    fill = color_wt,
    col = color_wt,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC2 KO 1
  atac_cdc2_ko_1 <- DataTrack(
    range = file.path(results_dir, "cDC2_KO_1.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC2_KO_1",
    fill = color_ko,
    col = color_ko,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC2 WT 2
  atac_cdc2_wt_2 <- DataTrack(
    range = file.path(results_dir, "cDC2_WT_2.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC2_WT_2",
    fill = color_wt,
    col = color_wt,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # cDC2 KO 2
  atac_cdc2_ko_2 <- DataTrack(
    range = file.path(results_dir, "cDC2_KO_2.bw"),
    genome = "mm10",
    type = "histogram",
    chromosome = chr,
    name = "cDC2_KO_2",
    fill = color_ko,
    col = color_ko,
    ylim = c(0, y_max),
    background.title = "white",
    col.axis = "black",
    col.title = "black",
    fontcolor.title = "black"
  )
  
  # Order: WT/KO pairs for each replicate
  track_list <- c(track_list, list(
    atac_cdc1_wt_1,
    atac_cdc1_ko_1,
    atac_cdc1_wt_2,
    atac_cdc1_ko_2,
    atac_cdc2_wt_1,
    atac_cdc2_ko_1,
    atac_cdc2_wt_2,
    atac_cdc2_ko_2
  ))
  
  # 5. Highlighted regions (differential peaks)
  if (gene_name == "Maged1") {
    highlight_gr <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = 94538000, end = 94548000)
    )
    highlight_color <- color_highlight_gray
  } else if (gene_name == "SpiB") {
    highlight_gr <- GRanges(
      seqnames = chr,
      ranges = IRanges(
        start = c(44526500, 44531000),
        end = c(44528500, 44532000)
      )
    )
    highlight_color <- color_highlight_pink
  } else {
    highlight_gr <- GRanges()
  }
  
  # Create highlight track if we have regions
  if (length(highlight_gr) > 0) {
    highlight_track <- AnnotationTrack(
      highlight_gr,
      name = "Diff Peaks",
      fill = highlight_color,
      col = highlight_color,
      alpha = 0.3
    )
    track_list <- c(track_list, list(highlight_track))
  }
  
  # 6. Plot - PNG format
  cat("Generating plot...\n")
  png(file.path(output_dir, output_file), 
      width = 12, 
      height = 14, 
      units = "in", 
      res = 300)
  plotTracks(
    track_list,
    from = start,
    to = end,
    chromosome = chr,
    background.title = "white",
    col.axis = "black",
    cex.axis = 0.8,
    cex.title = 0.9
  )
  dev.off()
  
  cat(paste0("Plot saved to: ", file.path(output_dir, output_file), "\n\n"))
}

# =============================================================================
# Generate the plots
# =============================================================================

# Panel D: Maged1 - ALL samples now
create_browser_plot(
  chr = "chr9",
  start = 94523802,
  end = 94553816,
  gene_name = "Maged1",
  output_file = "Maged1_browser_track.png",
  include_wt = TRUE
)

# Panel F: SpiB - All samples (WT + KO), grouped by replicate
create_browser_plot(
  chr = "chr7",
  start = 44525000,
  end = 44533000,
  gene_name = "SpiB",
  output_file = "SpiB_browser_track.png",
  include_wt = TRUE
)

cat("All plots generated successfully!\n")