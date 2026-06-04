# Final Project: ATAC-Seq Nextflow Pipeline Analysis

## Source of Data

**Publication:** De Sá Fernandes C, Novoszel P, Gastaldi T, Krauß D, Lang M, Rica R, Kutschat AP, Holcmann M, Ellmeier W, Seruggia D, Strobl H, Sibilia M. The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity. Cell Rep. 2024 Jun 25;43(6):114308. doi: 10.1016/j.celrep.2024.114308. Epub 2024 Jun 2. PMID: 38829740.

**Data Source:** GEO accession GSE266584

## Synopsis

Dendritic cells (DCs) are specialized immune cells that bridge innate and adaptive immunity by detecting pathogen and tumor antigens, and activating naive T cells in lymph nodes to initiate immune responses. There are two types of DC cells that play a large role in immune responses: conventional DC1 (cDC1) and conventional DC2 (cDC2). The development of cDCs is orchestrated by transcription factors (TFs) such as IRF4, IRF8, and SPIB. However, the epigenetic mechanisms regulating chromatin accessibility and gene expression during DC development are not completely understood. In this study, the researchers focused on the role of histone deacetylase 1 (HDAC1), which is an enzyme that modulates chromatin structure and TF accessibility, and has been shown to play a crucial role in the development of immune cells.

The study used ATAC-seq on cDC1 and cDC2 cells from mice, with wild type (WT) and HDAC1 knock-out (KO) samples. ATAC-seq data was used to map chromatin accessibility using Tn5 transposase, identifying regulatory elements such as enhancers and promoters. Linking the changes in chromatin landscape, gene expression, and identifying differentially accessible regions gave the researchers a better understanding of the transcriptional and epigenetic factors governing DC specification in development and function.

## Running the Pipeline

### Prerequisites

-   Nextflow (version 21.04 or later)
-   Docker or Singularity for container support
-   R (version 4.0 or later) with DiffBind and DESeq2 packages installed

### Setup

Clone this repository and navigate to the project directory. Configure the `nextflow.config` file with the following required input files:

-   `sra` - CSV file with sample ID and SRR number (format: sample,srr)
-   `genome` - Primary genome assembly (GRCm38.primary_assembly.genome.fa)
-   `gtf` - Primary assembly GTF file from GENCODE (gencode.vM10.primary_assembly.annotation.gtf)
-   `ucsc_genes` - UCSC list of TSS genes (mm10_genes.bed)
-   `blacklist` - List of blacklist genes (mm10-blacklist.v2.bed)
-   `adapters` - FASTA file containing adapter sequences for ATAC-seq data (adapters.fa)

### Execution

#### Step 1: Run the main ATAC-seq pipeline

``` bash
nextflow run main.nf -profile singularity,cluster
```

**Data Preprocessing (QC, Adapter Trimming, Alignment)**

-   **SRA_DOWNLOAD** (sratoolkit v3.2.1): Download raw ATAC-seq data from NCBI's SRA using fasterq-dump
-   **FASTQC_RAW** (fastqc v0.12.1): Quality assessment of raw FASTQ files
-   **TRIM** (trimmomatic v0.39): Adapter trimming and quality filtering with parameters `ILLUMINACLIP:2:20:7`, `LEADING:3`, `TRAILING:3`, `SLIDINGWINDOW:4:15`, `MINLEN:20`
-   **FASTQC_TRIM** (fastqc v0.12.1): Quality assessment of trimmed FASTQ files
-   **BOWTIE2_BUILD** (bowtie2 v2.5.4): Build genome index
-   **BOWTIE2_ALIGN** (bowtie2 v2.5.4): Align reads in `--very-sensitive` mode
-   **SAMTOOLS_FLAGSTAT** (samtools v1.21): Generate alignment statistics
-   **MULTIQC** (multiqc v1.25): Aggregate QC metrics into single HTML report

**Peak Calling**

-   **SAMTOOLS_SORT** (samtools v1.21): Sort BAM files by genomic coordinates
-   **SAMTOOLS_VIEW** (samtools v1.21): Filter BAM files (`-q 30`, `-F 4`, `grep -v chrM`)
-   **SAMTOOLS_IDX** (samtools v1.21): Index filtered BAM files
-   **TAGDIR** (homer v4.11): Create tag directories using makeTagDirectory
-   **FINDPEAKS** (homer v4.11): Call peaks with `-style histone`
-   **POS2BED** (homer v4.11): Convert peaks to BED format
-   **BEDTOOLS_MERGE** (bedtools v2.32.1): Generate consensus peak set
-   **BEDTOOLS_INTERSECT** (bedtools v2.32.1): Quantify read counts in peaks

**Visualization**

-   **BAMCOVERAGE** (deeptools v3.5.5): Convert BAM to bigWig format
-   **COMPUTEMATRIX** (deeptools v3.5.5): Compute signal matrix at TSS (`--referencePoint TSS`, `--beforeRegionStartLength 2000`, `--afterRegionStartLength 2000`, `--binSize 10`, `--missingDataAsZero`)

**ATAC-Seq QC Metrics**

-   **FRIP_SCORE** (bedtools v2.32.1): Calculate Fraction of Reads in Peaks
-   **TSS_ENRICHMENT** (homer v4.11): Calculate TSS enrichment scores (`-size 4000`, `-hist 10`)
-   **SUMMARIZE_TSS** (homer v4.11): Aggregate TSS enrichment results

#### Step 2: Run differential accessibility analysis

``` bash
Rscript -e diffbind.Rmd
```

Performs differential accessibility analysis using DiffBind with DESeq2: - Creates consensus peaks (`minOverlap = 2`) - TMM normalization (`score = DBA_SCORE_TMM_MINUS_FULL`) - Statistical testing with DESeq2 (`method = DBA_DESEQ2`, `bBlacklist = FALSE`, `bGreylist = FALSE`) - Filters results at p-value \< 0.01 - Categorizes peaks as "gain" (positive fold change) or "loss" (negative fold change)

**Outputs:** - `results/diffbind/cdc1_diffbind_results.csv` - `results/diffbind/cdc2_diffbind_results.csv` - `results/diffbind/cdc1_gain_peaks.bed` and `cdc1_loss_peaks.bed` - `results/diffbind/cdc2_gain_peaks.bed` and `cdc2_loss_peaks.bed`

#### Step 3: Run peak annotation and motif discovery

``` bash
nextflow run peak_analysis.nf -profile singularity,cluster
```

**Peak Annotation and Motif Discovery**

-   **FIND_MOTIFS_GENOME** (homer v4.11): Discover enriched motifs (`-size 200`, `-mask`, `-len 8,10,12`)
-   **ANNOTATE** (homer v4.11): Annotate peaks with genomic features using `-gtf`

**Differential Accessibility Visualization**

-   **COMPUTEMATRIX_DIFFERENTIAL** (deeptools v3.5.5): Compute matrices for gain/loss peaks
-   **HEATMAP** (deeptools v3.5.5): Generate heatmaps with `--groupLabels "Gained"` and `"Lost"`

#### Step 4: Extract gene lists for pathway enrichment

``` bash
Rscript annotpeaks_wfc.R
```

Annotates differential peaks with gene information and extracts gene lists for pathway enrichment analysis.

**Outputs:** - `results/cDC1_genes_for_enrichr.txt` - `results/cDC2_genes_for_enrichr.txt` - `results/cDC1_gain_genes.txt` and `cDC1_loss_genes.txt` - `results/cDC2_gain_genes.txt` and `cDC2_loss_genes.txt`

## Project Deliverables

This pipeline produces the following deliverables:

### 1. MultiQC Report

Comprehensive quality control report aggregating metrics from all preprocessing steps.

**Location:** `results/multiqc_report.html`

### 2. ATAC-seq Quality Metrics

-   **FRiP Scores:** Fraction of Reads in Peaks for each sample
-   **TSS Enrichment Scores:** Signal enrichment at transcription start sites

**Locations:** - `results/all_frip_scores.txt` - `results/tss_enrichment/` (individual TSS enrichment files) - Individual sample metrics: `results/{sample}.frip.txt` and `results/{sample}.flagstat.txt`

### 3. Consensus Peak Sets

Merged peaks across all samples for differential analysis.

**Location:** `results/consensus_peaks.bed`

### 4. Differential Accessibility Results

Statistical analysis results identifying significantly different peaks between WT and KO conditions.

**Locations:** - `results/diffbind/cdc1_diffbind_results.csv` (505 significant peaks) - `results/diffbind/cdc2_diffbind_results.csv` (265 significant peaks) - `results/diffbind/cdc1_gain_peaks.bed` and `cdc1_loss_peaks.bed` - `results/diffbind/cdc2_gain_peaks.bed` and `cdc2_loss_peaks.bed`

### 5. Peak Annotations

Genomic context and gene associations for differential peaks.

**Locations:** - `results/annotation/cDC1_annotated_peaks.txt` - `results/annotation/cDC2_annotated_peaks.txt` - `results/cDC1_annotated_with_fc.txt` - `results/cDC2_annotated_with_fc.txt`

### 6. Gene Lists for Enrichment Analysis

Unique gene symbols associated with differential peaks, ready for pathway enrichment tools like Enrichr.

**Locations:** - `results/cDC1_genes_for_enrichr.txt` - `results/cDC2_genes_for_enrichr.txt` - `results/cDC1_gain_genes.txt` and `cDC1_loss_genes.txt` - `results/cDC2_gain_genes.txt` and `cDC2_loss_genes.txt`

### 7. Motif Enrichment Results

Transcription factor binding motifs enriched in differential peaks.

**Locations:** - `results/cDC1_motifs/` (Homer motif analysis output) - `results/cDC2_motifs/` (Homer motif analysis output) - `results/homer_cdc1.jpg` (motif enrichment summary figure) - `results/homer_cdc2.jpg` (motif enrichment summary figure)

### 8. Visualization Files

-   **BigWig coverage tracks** for genome browser visualization
-   **Differential accessibility heatmaps** showing gain vs. loss peaks
-   **TSS enrichment heatmaps**
-   **Genome browser tracks** for specific loci (Maged1, Spib)
-   **Pathway enrichment figures**

**Locations:** - `results/*.bw` (8 bigWig files, one per sample) - `results/cDC1_differential_heatmap.png` - `results/cDC2_differential_heatmap.png` - `results/cDC1_tss_heatmap.png` - `results/cDC2_tss_heatmap.png` - `results/Maged1_browser_track.png` and `Maged1_browser_track.pdf` - `results/SpiB_browser_track.png` and `SpiB_browser_track.pdf` - `results/enrichr_cdc1.jpg` - `results/enrichr_cdc2.jpg`

### 9. Individual Sample Files

Per-sample outputs for quality control and downstream analysis: - Trimmed FASTQ files: `results/{sample}_trimmed.fastq.gz` - Filtered BAM files: `results/{sample}.filtered.bam` and `.bai` - Peak calls: `results/{sample}_peaks.txt` and `.bed` - Tag directories: `results/{sample}_tags/` - QC reports: `results/{sample}_trimmed_fastqc.html` - Trimmomatic logs: `results/{sample}_trimmomatic.log`

## Citation

If you use this pipeline, please cite the original publication:

De Sá Fernandes C, Novoszel P, Gastaldi T, et al. The histone deacetylase HDAC1 controls dendritic cell development and anti-tumor immunity. Cell Rep. 2024;43(6):114308.
