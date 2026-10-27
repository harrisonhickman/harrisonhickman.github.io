# Bulk RNA‑seq Analysis Project

## Overview
This project involved performing a complete bulk RNA‑seq analysis workflow using R and Bioconductor packages. The aim was to identify differentially expressed genes between two conditions, evaluate data quality, and interpret biological significance through enrichment analysis. This work demonstrates structured analytical workflow design, QC awareness, and biological interpretation aligned with clinical bioinformatics practice.

## Objectives
- Assess data quality and perform initial QC
- Normalise count data
- Identify differentially expressed genes (DEGs)
- Visualise biological variation
- Conduct GO enrichment analysis

## Workflow Summary
1. **Data import and QC**  
#### Navigate to sample directory
cd "data/sample_X"

#### Download FASTQ files from SRA
prefetch SRRXXXXXXX
fasterq-dump SRRXXXXXXX

#### Verify metadata
vdb-dump SRRXXXXXXX --rows 1

#### Quality control
run_fastqc.bat

#### Quantification with Kallisto
kallisto quant `
  -i "reference/human_transcriptome.idx" `
  -o kallisto_output `
  --single `
  -l 200 `
  -s 30 `
  SRRXXXXXXX.fastq

2. **Normalisation**  
   Applied appropriate normalisation methods to correct for sequencing depth and composition bias.

3. **Differential Expression**  
   Used statistical modelling to identify significantly up‑ and down‑regulated genes.

4. **Visualisation**  
   Generated PCA plots, volcano plots, heatmaps, and MA plots to explore biological patterns.

5. **Functional Enrichment**  
   Performed GO and KEGG analysis to identify pathways associated with DEGs.

## Key Results
(pending)

- **PCA plot** showing sample clustering  
- **Volcano plot** highlighting significant DEGs  
- **Heatmap** of top variable genes  
- **Enrichment barplot** showing GO/KEGG pathways  

## Tools & Packages
- R  
- DESeq2 / edgeR  
- clusterProfiler  
- ggplot2  
- tidyverse  

## Skills Demonstrated
- Bioinformatics workflow design  
- Statistical analysis of gene expression  
- Data QC and interpretation  
- Visualisation of biological variation  
- Reproducible analysis using R  
