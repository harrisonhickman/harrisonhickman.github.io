<img width="4032" height="1541" alt="image" src="https://github.com/user-attachments/assets/f9ecf16a-6eb8-4281-a125-fa78891ca371" /># Bulk RNA‑seq Analysis Project

## Overview
This project involved performing a complete bulk RNA‑seq analysis workflow using R and Bioconductor packages. The aim was to identify differentially expressed genes between two conditions, evaluate data quality, and interpret biological significance through enrichment analysis. This work demonstrates structured analytical workflow design, QC awareness, and biological interpretation aligned with clinical bioinformatics practice.

## Objectives
- Assess data quality and perform initial QC
- Normalise count data
- Identify differentially expressed genes (DEGs)
- Visualise biological variation
- Conduct GO enrichment analysis

## Melanoma RNA-seq Analysis Pipeline Overview
### 1. Data acquisition
•	Downloaded RNA-seq FASTQ files from SRA.
•	Dataset contained: 
o	3 × primary melanocyte replicates
o	3 × SK-Mel-28 replicates
o	3 × SK-Mel-147 replicates
o	3 × UACC-62 replicates
### 2. Transcript quantification
•	Built a Kallisto index using the reference transcriptome.
•	Quantified each sample individually using Kallisto.
•	Generated abundance.tsv files for all 12 samples.
### 3. Import into R
•	Imported all Kallisto outputs using tximport.
•	Combined all 12 samples into a single expression matrix.
•	Assigned sample names: 
o	Melanocytes_1-3
o	SK-Mel-28_1-3
o	SK-Mel-147_1-3
o	UACC-62_1-3
### 4. Metadata creation
•	Created a metadata table describing the experimental groups.
•	Experimental factor:
Plain Text
condition
├── Melanocytes
├── SK-Mel-28
├── SK-Mel-147
└── UACC-62
### 5. Differential expression analysis
•	Created a DESeq2 dataset using:
R
design = ~ condition
•	Performed differential expression analysis with DESeq2.
### 6. Comparisons performed
•	SK-Mel-28 vs Melanocytes
•	SK-Mel-147 vs Melanocytes
•	UACC-62 vs Melanocytes
### 7. Significant features identified
Using:R
padj < 0.05
Results:
Comparison	Significant transcripts
SK-Mel-28 vs Melanocytes	18,183
SK-Mel-147 vs Melanocytes	19,297
UACC-62 vs Melanocytes	18,176
### 8. Annotation of significant transcripts
•	Extracted significant transcripts.
•	Removed RefSeq version numbers.
•	Mapped RefSeq transcript IDs to gene symbols and gene names using:
R
org.Hs.eg.db
Example:
Plain Text
NM_001242460 → GHR
NM_003759 → SLC4A4
NM_012120 → CD2AP
### 9. Identification of top differentially expressed genes
For each comparison:
Upregulated genes
•	Sorted by decreasing log2 fold change.
•	Extracted top 20 upregulated genes.
Downregulated genes
•	Sorted by increasing log2 fold change.
•	Extracted top 20 downregulated genes.
### 10. Planned downstream analyses
•	PCA
•	Heatmap generation
•	Volcano plots
•	Identification of genes shared across all melanoma cell lines
•	Functional enrichment (GO / KEGG / Reactome)

## Key Results
### Top 10 Upregulated Genes
[Uploading image.png…]()



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
