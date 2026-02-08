# B-ALL Transcriptome Analysis Pipeline

## Overview
This repository contains the computational workflow used to identify molecular drivers of relapse in B-cell Acute Lymphoblastic Leukemia (B-ALL). The pipeline processes RNA-Seq count data to perform differential expression analysis, batch effect correction, functional annotation, and clinical survival validation.

## Repository Contents

### 1. `01_Master_DGE_Pipeline.R`
**Core Differential Expression Workflow**
* **Input:** Raw counts and clinical metadata.
* **Method:** Uses `DESeq2` for normalization and statistical testing.
* **Output:** 
    * Volcano plots identifying significant up/down-regulated genes.
    * Annotated Heatmaps of top 50 differentially expressed genes.
    * CSV exports of all statistical results.

### 2. `02_Gradient_Pattern_Analysis.R`
**Longitudinal/Gradient Discovery**
* **Goal:** Identifies genes that follow a specific "staircase" expression pattern across three clinical states (e.g., Favourable → Intermediate → Unfavourable).
* **Method:** Uses Likelihood Ratio Testing (LRT) and custom boolean logic to classify genes into "Staircase Up", "Staircase Down", or "Late-Stage Specific" patterns.

### 3. `03_Batch_Correction_QC.R`
**Quality Control & Batch Effect Removal**
* **Method:** Implements `ComBat-seq` (sva package) to remove technical batch effects while preserving biological variation.
* **Validation:** Generates PCA and UMAP plots before and after correction to visualize data integrity.

### 4. `04_Biotype_Annotation.R`
**Functional Annotation & ncRNA Discovery**
* **Goal:** Classifies differentially expressed genes into biological categories (Protein-coding, lncRNA, miRNA).
* **Method:** Queries Ensembl BioMart to retrieve gene biotypes.
* **Visualization:** Automatically generates bar plots for the top 10 up- and down-regulated long non-coding RNAs (lncRNAs) and miRNAs.

### 5. `05_Survival_Analysis_KaplanMeier.R`
**Clinical Validation**
* **Goal:** Correlates gene expression levels with patient survival outcomes.
* **Method:** Stratifies patients into "High" vs "Low" expression groups (median split) and performs Kaplan-Meier analysis with Log-Rank testing.
* **Output:** Publication-ready survival curves with risk tables and p-values.

## Dependencies
The workflow relies on R (v4.0+) and the following key packages:
* **Bioconductor:** `DESeq2`, `ComplexHeatmap`, `sva`, `org.Hs.eg.db`, `biomaRt`
* **Visualization:** `ggplot2`, `pheatmap`, `ggrepel`, `ggpubr`, `survminer`
* **Statistics & Data:** `survival`, `dplyr`, `tidyr`, `tibble`

## Author
**Gadha K Leons** *PhD Candidate, Cancer Genomics & Molecular Biology* *All India Institute of Medical Sciences (AIIMS), New Delhi*