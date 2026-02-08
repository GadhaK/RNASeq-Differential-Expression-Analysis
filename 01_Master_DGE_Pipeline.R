# MASTER DGE PIPELINE: Differential Expression, Volcano, Heatmap & Pathways


# 1. SETUP

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(DESeq2, ggplot2, pheatmap, RColorBrewer, ggrepel, clusterProfiler, 
               enrichplot, org.Hs.eg.db, dplyr, biomaRt, ggvenn)

# 2. LOAD DATA

counts_file <- "merged_counts_BALL.csv"
metadata_file <- "merged_metadata_BALL.csv"

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1)
metadata <- metadata[colnames(counts), ] # Align samples

# 3. DESEQ2 OBJECT

# NOTE: If correcting for batch, use: design = ~ Batch + Relapse_status
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Relapse_status)
dds$Relapse_status <- relevel(dds$Relapse_status, ref = "Baseline") # Set reference

# Pre-filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run Analysis
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE) # For plotting

# 4. AUTOMATED ANALYSIS FUNCTION
# ----------------------------------------------------------------------------
run_comparison <- function(contrast, prefix) {

  # A. Get Results
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  res$symbol[is.na(res$symbol)] <- rownames(res)[is.na(res$symbol)]
  
  # B. Save CSV
  res_df <- as.data.frame(res) %>% arrange(padj)
  write.csv(res_df, paste0(prefix, "_All_Results.csv"))
  
  # C. Volcano Plot
  res_df$reg <- ifelse(res_df$log2FoldChange > 1.5 & res_df$padj < 0.05, "UP",
                       ifelse(res_df$log2FoldChange < -1.5 & res_df$padj < 0.05, "DOWN", "NS"))
  
  p_vol <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=reg)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("UP"="firebrick3", "DOWN"="navy", "NS"="grey")) +
    geom_text_repel(data=head(subset(res_df, reg!="NS"), 20), aes(label=symbol)) +
    theme_minimal() + ggtitle(paste("Volcano:", prefix))
  ggsave(paste0(prefix, "_Volcano.png"), p_vol)
  
  # D. Heatmap (Top 50)
  top50 <- head(subset(res_df, reg!="NS")$symbol, 50)
  if(length(top50) > 0) {
    # Match symbols back to Ensembl IDs for the matrix
    top_ids <- rownames(res_df[res_df$symbol %in% top50, ])[1:50] 
    mat <- assay(vsd)[top_ids, ]
    rownames(mat) <- res_df[top_ids, "symbol"]
    
    # Create Annotation DataFrame
    annotation_col <- as.data.frame(colData(dds)[, c("Relapse_status")])
    rownames(annotation_col) <- colnames(dds)
    colnames(annotation_col) <- "Group" # Rename for the legend title
    
    # Generate Heatmap
    pheatmap(mat, 
             scale = "row", 
             cluster_cols = TRUE, 
             annotation_col = annotation_col,
             show_colnames = FALSE,
             main = paste("Top 50:", prefix), 
             filename = paste0(prefix, "_Heatmap.png"))
  }
}

# 5. RUN COMPARISONS

run_comparison(c("Relapse_status", "Relapse", "Baseline"), "Relapse_vs_Baseline")

# 6. CANDIDATE GENE EXPRESSION (Boxplot)

library(ggpubr)
library(ggplot2)

plot_gene <- function(gene_name) {
  
  # 1. Map ID (Same as before)
  gene_id <- mapIds(org.Hs.eg.db, keys = gene_name, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  
  if(is.na(gene_id)) return(message("Gene not found"))
  if(!gene_id %in% rownames(dds)) return(message("Gene filtered out"))
  
  # 2. Get Data
  d <- plotCounts(dds, gene=gene_id, intgroup="Relapse_status", returnData=TRUE)
  
  # 3. Plot
  p <- ggplot(d, aes(x=Relapse_status, y=count, fill=Relapse_status)) +
    geom_boxplot(alpha=0.6, outlier.shape=NA) + 
    geom_jitter(width=0.2, size=1.5) +
    scale_y_log10() + 
    labs(title = paste(gene_name, "Expression"), 
         y = "Normalized Counts (Log10)") +
    theme_bw() +
    theme(legend.position = "none")
  
  # 4. Add Statistics Automatically
  # If you have 2 groups, this adds a T-test.
  # If you have >2 groups, this adds a Kruskal-Wallis global p-value.
  p <- p + stat_compare_means(method = "t.test") 
  
  print(p)
}

# Usage Examples:
plot_gene("CRLF2")
plot_gene("IKZF1")
