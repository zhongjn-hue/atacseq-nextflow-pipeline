ATAC-seq Analysis of Chromatin Accessibility in cDC1 and cDC2
1. Data Source and Publication

This project is based on publicly available ATAC-seq data generated to study chromatin accessibility changes in conventional dendritic cell subsets (cDC1 and cDC2) under wild-type (WT) and knockout (KO) conditions.

The original dataset and biological context are described in the following publication:
De Sá Fernandes, C., Novoszel, P., Gastaldi, T., Krauß, D., Lang, M., Rica, R., Kutschat, A. P., Holcmann, M., Ellmeier, W., Seruggia, D., Strobl, H., & Sibilia, M. (2024). The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity. Cell reports, 43(6), 114308. https://doi.org/10.1016/j.celrep.2024.114308


2. Repository Contents

This repository contains the entire analysis workflow and code, including:

-A Nextflow pipeline for ATAC-seq preprocessing
-An R Markdown report for downstream analysis, visualization, and figure reproduction
-Configuration files and module definitions required to reproduce the workflow

No raw FASTQ files, BAM files, or large output files are included, in accordance with course requirements.

3. ATAC-seq Processing Pipeline

ATAC-seq preprocessing is implemented using Nextflow, with modular steps for each stage of the analysis:

-Quality control of raw reads (FastQC)
-Adapter trimming and quality filtering (Trimmomatic)
-Alignment to the mouse reference genome (mm10) using Bowtie2
-Removal of mitochondrial reads
-BAM sorting and indexing (samtools)
-Peak calling (MACS2 / BEDTools)
-Generation of normalized coverage tracks (bigWig)
-Aggregation of quality control metrics (MultiQC)

The workflow is designed to be fully reproducible and portable.

4. How to Run the Pipeline
Requirements:
-Nextflow
-Java
-Docker / Singularity or a compatible HPC environment
-R (v4.4.3) with required Bioconductor packages

Running the Nextflow Pipeline
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --outdir results


All parameters related to genome reference files and computational resources are specified in the nextflow.config file.

Running Downstream Analysis

Downstream statistical analysis and visualization are performed using an R Markdown document:

analysis.Rmd


This file contains all R code used for:

Differential chromatin accessibility analysis (DiffBind / DESeq2)

Peak annotation (ChIPseeker)

Motif enrichment analysis (chromVAR, motifmatchr)

TSS-centered ATAC-seq heatmap visualization

Reproduction of ATAC-seq results corresponding to Figure 6A–B of the original study

The R Markdown can be rendered using:

rmarkdown::render("analysis.Rmd")

5. Project Deliverables

The following deliverables are included in this repository, as required for the final project:

✔ A fully reproducible ATAC-seq preprocessing pipeline implemented in Nextflow

✔ An R Markdown analysis notebook containing all downstream analysis code

✔ ATAC-seq–specific quality control metrics (TSS-centered enrichment heatmaps)

✔ Reproduction of key ATAC-seq accessibility patterns from the original publication

✔ A complete README documenting data sources, pipeline execution, and deliverables

All analyses are fully scripted and reproducible without inclusion of raw data files.

6. Software Versions

Key tools and packages used in this project include:

-FastQC v0.12.1

-Trimmomatic v0.39

-Bowtie2 v2.5.4

-samtools v1.21

-MACS2 v2.1.2

-BEDTools v2.31.1

-MultiQC v1.25

-R v4.4.3

-DiffBind, ChIPseeker, chromVAR, motifmatchr, DESeq2