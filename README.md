# B-ALL Transcriptome Analysis Pipeline

## Overview
This repository contains the computational workflow used to identify molecular drivers of relapse in B-cell Acute Lymphoblastic Leukemia (B-ALL). The pipeline processes RNA-Seq count data to perform differential expression analysis, batch effect correction, and pattern recognition for genes associated with disease progression.

## Repository Contents

### 1. `01_Master_DGE_Pipeline.R`
**Core Differential Expression Workflow**
  * **Input:** Raw counts and clinical metadata.
* **Method:** Uses `DESeq2` for normalization and statistical testing.
* **Output:** * Volcano plots identifying significant up/down-regulated genes.
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

## Dependencies
* **Bioconductor:** `DESeq2`, `ComplexHeatmap`, `sva`, `org.Hs.eg.db`
* **Visualization:** `ggplot2`, `pheatmap`, `ggrepel`, `ggpubr`
* **Data Wrangling:** `dplyr`, `tidyr`, `tibble`

## Author
**Gadha K Leons** *PhD Candidate, Cancer Genetics & Molecular Biology* *All India Institute of Medical Sciences (AIIMS), New Delhi*