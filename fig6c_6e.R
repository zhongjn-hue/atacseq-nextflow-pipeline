# integrate_atac_rna_plot.R
library(tidyverse)
library(ggrepel)

# Read ATAC data (already annotated with fold changes)
cdc1_atac <- read.delim("results/cDC1_annotated_with_fc.txt", header = TRUE)
cdc2_atac <- read.delim("results/cDC2_annotated_with_fc.txt", header = TRUE)

# Read RNA-seq data
cdc1_rna <- read.delim("data/rna_seq/cDC1_log2FC.tsv", header = TRUE)
cdc2_rna <- read.delim("data/rna_seq/cDC2_log2FC.tsv", header = TRUE)

# Genes to highlight (from the paper)
cdc1_genes <- c("Nectin2", "Il13ra1", "Maged1", "Rasa1", "Nfil3")
cdc2_genes <- c("Spib", "Irf8", "Susd2", "Nectin2", "Csf3r")

# Function to integrate ATAC and RNA
integrate_atac_rna <- function(atac_df, rna_df, genes_to_label, cell_type) {
  
  # Summarize ATAC by gene (average if multiple peaks per gene)
  atac_summary <- atac_df %>%
    filter(!is.na(gene_name) & gene_name != "") %>%
    group_by(gene_name) %>%
    summarize(
      ATAC_log2FC = mean(log2FoldChange, na.rm = TRUE),
      ATAC_padj = min(padj, na.rm = TRUE),
      n_peaks = n()
    ) %>%
    ungroup()
  
  # Merge ATAC and RNA by gene name
  integrated <- inner_join(
    atac_summary,
    rna_df,
    by = c("gene_name" = "gene_id")
  )
  
  # Color based on which quadrant the point falls in
  integrated <- integrated %>%
    mutate(
      quadrant = case_when(
        ATAC_log2FC > 0 & log2FC > 0 ~ "upper_right",  # RED
        ATAC_log2FC < 0 & log2FC < 0 ~ "lower_left",   # BLUE
        TRUE ~ "other"  # GRAY
      ),
      highlight = gene_name %in% genes_to_label
    )
  
  cat("Cell type:", cell_type, "\n")
  cat("Total genes matched:", nrow(integrated), "\n")
  cat("Upper right (red):", sum(integrated$quadrant == "upper_right"), "\n")
  cat("Lower left (blue):", sum(integrated$quadrant == "lower_left"), "\n")
  cat("Other (gray):", sum(integrated$quadrant == "other"), "\n")
  cat("Highlighted genes found:", sum(integrated$highlight), "\n\n")
  
  return(integrated)
}

# Integrate data
cdc1_integrated <- integrate_atac_rna(cdc1_atac, cdc1_rna, cdc1_genes, "cDC1")
cdc2_integrated <- integrate_atac_rna(cdc2_atac, cdc2_rna, cdc2_genes, "cDC2")

# Plot cDC1
png("results/cDC1_ATAC_RNA_integration.png", width = 8, height = 8, units = "in", res = 300)
ggplot(cdc1_integrated, aes(x = log2FC, y = ATAC_log2FC)) +
  # Add quadrant lines (thick black lines at x=0 and y=0)
  geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  # Plot points colored by quadrant
  geom_point(aes(color = quadrant), alpha = 0.6, size = 2) +
  geom_point(data = filter(cdc1_integrated, highlight),
             aes(color = quadrant), size = 3, shape = 21, stroke = 1.5) +
  # Add gene labels
  geom_text_repel(
    data = filter(cdc1_integrated, highlight),
    aes(label = gene_name),
    size = 4,
    fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c("upper_right" = "red", "lower_left" = "blue", "other" = "gray70"),
    guide = "none"
  ) +
  labs(
    title = "cDC1",
    x = "log2FoldChange RNA",
    y = "log2FoldChange ATAC"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
dev.off()

# Plot cDC2
png("results/cDC2_ATAC_RNA_integration.png", width = 8, height = 8, units = "in", res = 300)
ggplot(cdc2_integrated, aes(x = log2FC, y = ATAC_log2FC)) +
  # Add quadrant lines (thick black lines at x=0 and y=0)
  geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  # Plot points colored by quadrant
  geom_point(aes(color = quadrant), alpha = 0.6, size = 2) +
  geom_point(data = filter(cdc2_integrated, highlight),
             aes(color = quadrant), size = 3, shape = 21, stroke = 1.5) +
  # Add gene labels
  geom_text_repel(
    data = filter(cdc2_integrated, highlight),
    aes(label = gene_name),
    size = 4,
    fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c("upper_right" = "red", "lower_left" = "blue", "other" = "gray70"),
    guide = "none"
  ) +
  labs(
    title = "cDC2",
    x = "log2FoldChange RNA",
    y = "log2FoldChange ATAC"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
dev.off()

# Print highlighted genes for verification
cat("=== cDC1 Highlighted Genes ===\n")
print(cdc1_integrated %>% 
        filter(highlight) %>% 
        dplyr::select(gene_name, ATAC_log2FC, log2FC, quadrant))

cat("\n=== cDC2 Highlighted Genes ===\n")
print(cdc2_integrated %>% 
        filter(highlight) %>% 
        dplyr::select(gene_name, ATAC_log2FC, log2FC, quadrant))

cat("\nPlots saved to results/\n")