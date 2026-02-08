# BATCH CORRECTION & QC PIPELINE

library(sva)
library(DESeq2)
library(ggplot2)
library(umap)

# 1. LOAD DATA
counts <- read.delim("raw_counts.txt", row.names=1, check.names=FALSE)
meta <- read.delim("metadata.txt", row.names=1)

# Filter low expression
counts <- counts[rowSums(counts >= 10) >= ncol(counts)*0.9, ]

# 2. BATCH CORRECTION (ComBat-seq)

# Use this if you have raw counts and want to preserve integer nature for DESeq2
adjusted_counts <- ComBat_seq(as.matrix(counts), 
                              batch = meta$Batch, 
                              group = meta$Risk_status) # Protect the biological variable!

write.csv(adjusted_counts, "combat_seq_corrected_counts.csv")

# 3. VISUALIZATION (PCA & UMAP)
# Compare Before vs After

plot_pca_custom <- function(cnts, metadata, title) {
  dds_temp <- DESeqDataSetFromMatrix(cnts, metadata, ~1)
  vsd_temp <- vst(dds_temp, blind=TRUE)
  pcaData <- plotPCA(vsd_temp, intgroup=c("Risk_status", "Batch"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=Risk_status, shape=Batch)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    ggtitle(title) + theme_bw()
}

# Plot Original
p1 <- plot_pca_custom(counts, meta, "Before Correction")
ggsave("PCA_Before_Batch.png", p1)

# Plot Corrected
p2 <- plot_pca_custom(adjusted_counts, meta, "After ComBat-seq")
ggsave("PCA_After_Batch.png", p2)