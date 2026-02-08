# BIOTYPE ANNOTATION & VISUALIZATION
# Annotate DEGs with biotypes (coding vs. non-coding) and visualize top ncRNAs.

# SETUP
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(biomaRt, dplyr, ggplot2, tidyverse)

input_file <- "Relapse_vs_Baseline_All_Results.csv" 

# LOAD DATA
degs <- read.csv(input_file, header = TRUE)
colnames(degs)[1] <- "ensembl_gene_id" # Ensure ID column matches

# QUERY BIOMART
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                   filters = "ensembl_gene_id",
                   values = degs$ensembl_gene_id,
                   mart = mart)

# MERGE & CLASSIFY
annotated_degs <- left_join(degs, gene_info, by = "ensembl_gene_id")

# Replace empty gene names with Ensembl IDs
annotated_degs$external_gene_name <- ifelse(annotated_degs$external_gene_name == "", 
                                            annotated_degs$ensembl_gene_id, 
                                            annotated_degs$external_gene_name)

# Split by Category
mRNA   <- annotated_degs %>% filter(gene_biotype == "protein_coding")
lncRNA <- annotated_degs %>% filter(gene_biotype == "lncRNA")
miRNA  <- annotated_degs %>% filter(gene_biotype == "miRNA")

# 6. EXPORT TABLES
write.csv(mRNA,   "DE_Protein_Coding_Genes.csv", row.names = FALSE)
write.csv(lncRNA, "DE_lncRNAs.csv", row.names = FALSE)
write.csv(miRNA,  "DE_miRNAs.csv", row.names = FALSE)

# VISUALIZATION (Top Up/Down Bar Plot)

plot_top_biotype <- function(df, biotype_name, top_n = 10) {
  
  if(nrow(df) < 2) return(message(paste("Not enough genes to plot", biotype_name)))
  
  # Select Top Up and Top Down
  top_up <- df %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) %>% head(top_n)
  top_down <- df %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange) %>% head(top_n)
  
  plot_data <- bind_rows(top_up, top_down)
  
  # Create "Regulation" column for coloring
  plot_data$Regulation <- ifelse(plot_data$log2FoldChange > 0, "Upregulated", "Downregulated")
  
  # Create Plot
  p <- ggplot(plot_data, aes(x = reorder(external_gene_name, log2FoldChange), 
                             y = log2FoldChange, 
                             fill = Regulation)) +
    geom_col() +
    coord_flip() +  # Flip to make bars horizontal
    scale_fill_manual(values = c("Upregulated" = "firebrick3", "Downregulated" = "navy")) +
    labs(title = paste("Top DE", biotype_name),
         x = "Gene Symbol",
         y = "Log2 Fold Change") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
  
  # Save Plot
  filename <- paste0("Top_", biotype_name, "_Barplot.png")
  ggsave(filename, p, width = 7, height = 6)
  print(p)
}

# Run the function for lncRNA and miRNA
plot_top_biotype(lncRNA, "lncRNAs")
plot_top_biotype(miRNA, "miRNAs")
