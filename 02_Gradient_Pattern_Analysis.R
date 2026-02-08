# GRADIENT PATTERN ANALYSIS (The "Staircase" Effect)

# Focus: Identifying genes that show a progressive change across 3 groups
# (e.g., Favourable -> Intermediate -> Unfavourable OR Favourable -> Int -> Unfav)

# 1. SETUP & LIBRARIES
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(DESeq2, ggplot2, dplyr, tidyr, tibble, pheatmap, org.Hs.eg.db, RColorBrewer)


counts_file   <- "merged_counts_BALL.csv"
metadata_file <- "merged_metadata_BALL.csv"

group_col <- "Prognosis" 

# Define the logical order of your groups (Low -> Med -> High severity)
group_levels <- c("Favourable", "Intermediate", "Unfavourable") 

# 3. DATA LOADING & DESEQ2
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1)

common <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common]
metadata <- metadata[common, ]

# Create DESeq Object with LRT
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = as.formula(paste("~", group_col)))

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set Factor Levels (Critical for "Staircase" logic)
colData(dds)[[group_col]] <- factor(colData(dds)[[group_col]], levels = group_levels)

# Run DESeq2 Likelihood Ratio Test (LRT)
# This finds genes that change anywhere across the groups

dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
res <- results(dds)

# 4. PATTERN CLASSIFICATION LOGIC

# A. Extract Significant Genes (padj < 0.05)
sig_genes_df <- as.data.frame(res) %>% 
  rownames_to_column("Gene") %>% 
  filter(padj < 0.05) %>% 
  arrange(padj)

# Map IDs to Symbols
sig_genes_df$Symbol <- mapIds(org.Hs.eg.db, keys=sig_genes_df$Gene, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
sig_genes_df$Symbol[is.na(sig_genes_df$Symbol)] <- sig_genes_df$Gene[is.na(sig_genes_df$Symbol)]

# B. Calculate Mean Expression per Group
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)

# Subset to significant genes
sig_counts_long <- norm_counts[sig_genes_df$Gene, ] %>%
  as.data.frame() %>% 
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "SampleID", values_to = "Expression") %>%
  left_join(rownames_to_column(as.data.frame(colData(dds)), "SampleID"), by = "SampleID")

# Calculate Group Means
gene_means <- sig_counts_long %>%
  group_by(Gene, .data[[group_col]]) %>%
  summarise(Mean = mean(Expression), .groups = 'drop') %>%
  pivot_wider(names_from = all_of(group_col), values_from = Mean)

# C. Define Patterns (Using the 3 levels defined in CONFIG)
A <- group_levels[1] # e.g. Favourable
B <- group_levels[2] # e.g. Intermediate
C <- group_levels[3] # e.g. Unfavourable

classified_genes <- gene_means %>%
  mutate(
    # Staircase UP: A < B < C
    Pattern_Up = (.data[[B]] > .data[[A]]) & (.data[[C]] > .data[[B]]),
    
    # Staircase DOWN: A > B > C
    Pattern_Down = (.data[[A]] > .data[[B]]) & (.data[[B]] > .data[[C]]),
    
    # Specific to Last Stage: (C is highest, A and B are lower)
    Pattern_LateHigh = (.data[[C]] > .data[[B]]) & (.data[[C]] > .data[[A]])
  ) %>%
  left_join(sig_genes_df[, c("Gene", "Symbol")], by = "Gene") # Add symbols back

# D. Select Candidates
final_candidates <- classified_genes %>% filter(Pattern_Up == TRUE)

if(nrow(final_candidates) == 0) {
  message("No perfect 'Staircase Up' genes found. Switching to 'Late High' pattern.")
  final_candidates <- classified_genes %>% filter(Pattern_LateHigh == TRUE)
}

message(paste("Found", nrow(final_candidates), "candidate genes matching the pattern."))

# 5. VISUALIZATION

# Boxplot of Top 4 Candidates
top_genes <- head(final_candidates$Gene, 4)
plot_data <- sig_counts_long %>% filter(Gene %in% top_genes)

# Add Symbol for Plot Label
plot_data <- plot_data %>% left_join(sig_genes_df[,c("Gene", "Symbol")], by="Gene")

p1 <- ggplot(plot_data, aes(x = .data[[group_col]], y = Expression, fill = .data[[group_col]])) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Symbol, scales = "free_y") +
  theme_bw() +
  labs(title = "Top Gradient Candidates", y = "VST Normalized Expression") +
  theme(legend.position = "none")

print(p1)
ggsave("Gradient_Candidates_Boxplot.png", p1, width = 8, height = 6)


# Heatmap of All Pattern-Matching Genes
# Prepare Matrix
heatmap_genes <- final_candidates$Gene
if(length(heatmap_genes) > 50) heatmap_genes <- head(heatmap_genes, 50) # Cap at 50 for readability

mat_scaled <- norm_counts[heatmap_genes, ]
# Scale rows (Z-score)
mat_scaled <- t(scale(t(mat_scaled)))

# Replace rownames with Symbols
rownames(mat_scaled) <- final_candidates$Symbol[match(rownames(mat_scaled), final_candidates$Gene)]

# Annotation Bar
anno_col <- as.data.frame(colData(dds)[, group_col, drop=FALSE])
colnames(anno_col) <- "Group"

# Sort matrix columns by group
ordered_samples <- rownames(anno_col)[order(anno_col$Group)]
mat_scaled <- mat_scaled[, ordered_samples]
anno_col <- anno_col[ordered_samples, , drop=FALSE]

# 3. Plot
pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE, # Keep our sorted order
         annotation_col = anno_col,
         show_colnames = FALSE,
         scale = "none", # Already scaled
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Gradient Pattern Heatmap",
         filename = "Gradient_Pattern_Heatmap.png")

# 6. EXPORT RESULTS
write.csv(classified_genes, "Pattern_Classification_Analysis.csv", row.names = FALSE)
write.csv(final_candidates, "Selected_Gradient_Candidates.csv", row.names = FALSE)