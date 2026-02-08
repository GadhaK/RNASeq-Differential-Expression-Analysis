# SURVIVAL ANALYSIS (Kaplan-Meier)
# Correlate gene expression levels (High vs Low) with patient survival.

# SETUP
pacman::p_load(survival, survminer, ggplot2, dplyr)

counts_file <- "merged_counts_length.csv" # Needs a 'Length' column for TPM
meta_file   <- "merged_metadata_BALL.csv"

target_gene_id   <- "ENSG00000066336" 
target_gene_name <- "SPI1" 

# Clinical columns in metadata
time_col   <- "EF_duration"                # Time in days/months
status_col <- "Event_0lost1relapsedeath"   # 1=Event, 0=Censored

# DATA PROCESSING (Counts -> TPM)
counts_l <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta     <- read.csv(meta_file, row.names = 1)

# Calculate TPM (Transcripts Per Million)
gene_lengths <- counts_l$Length
counts_only  <- counts_l[, -ncol(counts_l)]

rpk <- counts_only / (gene_lengths / 1000)
tpm <- t(t(rpk) / (colSums(rpk) / 1e6))
log_tpm <- log2(tpm + 1)

# Align Metadata and Expression
common <- intersect(colnames(log_tpm), rownames(meta))
log_tpm <- log_tpm[, common]
meta    <- meta[common, ]

# PREPARE SURVIVAL GROUPS
# Extract expression for the target gene
if(!target_gene_id %in% rownames(log_tpm)) stop("Gene ID not found in matrix!")

gene_values <- as.numeric(log_tpm[target_gene_id, ])
cutoff      <- median(gene_values) # Median split

meta$Expression_Group <- ifelse(gene_values >= cutoff, "High", "Low")
meta$Expression_Group <- factor(meta$Expression_Group, levels = c("Low", "High"))

# 5. RUN SURVIVAL MODEL
surv_obj <- Surv(time = meta[[time_col]], event = meta[[status_col]])
fit      <- survfit(surv_obj ~ Expression_Group, data = meta)

# 6. PLOT
p <- ggsurvplot(
  fit,
  data = meta,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  title = paste("Survival Analysis:", target_gene_name),
  xlab = "Time (Days)",
  legend.labs = c("Low", "High"),
  ggtheme = theme_classic()
)

print(p)

# Save
ggsave(paste0("Survival_", target_gene_name, ".png"), print(p), width = 6, height = 5)
