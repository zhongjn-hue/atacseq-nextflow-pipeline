# integrate_atac_with_genes.R
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(tidyverse)  # Load tidyverse AFTER Bioconductor packages

# Read DiffBind results
cat("Reading DiffBind results...\n")
cdc1_diffbind <- read.csv("results/diffbind/cdc1_diffbind_results.csv") %>%
  filter(significant == TRUE)

cdc2_diffbind <- read.csv("results/diffbind/cdc2_diffbind_results.csv") %>%
  filter(significant == TRUE)

cat("cDC1 significant peaks:", nrow(cdc1_diffbind), "\n")
cat("cDC2 significant peaks:", nrow(cdc2_diffbind), "\n")

# Function to annotate and format
annotate_peaks <- function(diffbind_df, cell_type) {
  cat("\nAnnotating", cell_type, "peaks...\n")
  
  # Convert to GRanges
  gr <- makeGRangesFromDataFrame(diffbind_df, 
                                 seqnames.field = "seqnames",
                                 start.field = "start", 
                                 end.field = "end",
                                 keep.extra.columns = TRUE)
  
  # Annotate with ChIPseeker
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  annotated <- annotatePeak(gr, 
                            tssRegion = c(-3000, 3000), 
                            TxDb = txdb, 
                            annoDb = "org.Mm.eg.db")
  
  # Convert to dataframe
  anno_df <- as.data.frame(annotated)
  
  # Format output using dplyr::select explicitly
  result <- anno_df %>%
    dplyr::select(seqnames, start, end, width, strand,
                  Conc, Conc_KO, Conc_WT, Fold, p.value, FDR,
                  annotation, geneChr, geneStart, geneEnd, geneLength, 
                  geneStrand, geneId, transcriptId, distanceToTSS, SYMBOL, GENENAME) %>%
    dplyr::rename(gene_name = SYMBOL,
                  gene_description = GENENAME,
                  log2FoldChange = Fold,
                  padj = FDR)
  
  cat("Annotated", nrow(result), "peaks with", 
      sum(!is.na(result$gene_name)), "gene assignments\n")
  
  return(result)
}

# Annotate both cell types
cdc1_annotated <- annotate_peaks(cdc1_diffbind, "cDC1")
cdc2_annotated <- annotate_peaks(cdc2_diffbind, "cDC2")

# Save results
cat("\nSaving results...\n")
write.table(cdc1_annotated, "results/cDC1_annotated_with_fc.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(cdc2_annotated, "results/cDC2_annotated_with_fc.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done! Files saved to results/\n")

# Preview
cat("\n=== cDC1 Preview ===\n")
print(head(cdc1_annotated %>% dplyr::select(seqnames, start, end, gene_name, annotation, log2FoldChange, padj)))

cat("\n=== cDC2 Preview ===\n")
print(head(cdc2_annotated %>% dplyr::select(seqnames, start, end, gene_name, annotation, log2FoldChange, padj)))

# Extract unique gene names for Enrichr
cat("\n=== Extracting gene lists for Enrichr ===\n")

# cDC1 genes
cdc1_genes <- cdc1_annotated %>%
  filter(!is.na(gene_name)) %>%
  pull(gene_name) %>%
  unique() %>%
  sort()

cat("cDC1 unique genes:", length(cdc1_genes), "\n")

# cDC2 genes
cdc2_genes <- cdc2_annotated %>%
  filter(!is.na(gene_name)) %>%
  pull(gene_name) %>%
  unique() %>%
  sort()

cat("cDC2 unique genes:", length(cdc2_genes), "\n")

# Save gene lists (one gene per line for Enrichr)
writeLines(cdc1_genes, "results/cDC1_genes_for_enrichr.txt")
writeLines(cdc2_genes, "results/cDC2_genes_for_enrichr.txt")

cat("\nGene lists saved:\n")
cat("  - results/cDC1_genes_for_enrichr.txt\n")
cat("  - results/cDC2_genes_for_enrichr.txt\n")
