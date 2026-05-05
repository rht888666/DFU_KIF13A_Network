# ============================================================================
# DFU External Validation - Multi-Gene Validation (High-Impact Style)
# Script 33_ExtVal_MultiGene: Generate Advanced Figure 7A
# Target: Validate KIF13A + Its Team (Network Module)
# ============================================================================

rm(list = ls())
gc()

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
validation_dir <- file.path(results_dir, "Validation_Data")
gse_dir <- file.path(validation_dir, "GSE134431")
xlsx_files <- list.files(gse_dir, pattern = "\\.xlsx$", full.names = TRUE)
target_file <- xlsx_files[1]

# 1. Robust Header Detection (Reusing Logic)
final_df <- NULL
# Function to check cols
check_cols <- function(df) {
    cols <- colnames(df)
    has_logfc <- any(grepl("Log2|Fold|FC|logFC", cols, ignore.case = TRUE))
    has_pval <- any(grepl("P.Value|PValue|p-val|adj.P.Val", cols, ignore.case = TRUE))
    return(has_logfc && has_pval)
}
# Try skips
for (skip_n in 0:10) {
    try(
        {
            df <- read_excel(target_file, skip = skip_n, n_max = 5)
            if (check_cols(df)) {
                final_df <- read_excel(target_file, skip = skip_n)
                break
            }
        },
        silent = TRUE
    )
}
if (is.null(final_df)) stop("Header detection failed.")

# 2. Extract Data
cols <- colnames(final_df)
logfc_col <- cols[grep("Log2|Fold|FC|logFC", cols, ignore.case = TRUE)][1]
pval_col <- cols[grep("P.Value|PValue|p-val|adj.P.Val", cols, ignore.case = TRUE)][1]
gene_col <- cols[1]

# 3. Define the "Transport Module"
target_genes <- c("KIF13A", "EPN1", "CLIP1", "RAB11A", "RAB4A", "CTNNB1", "CDH1")

# Clean Gene Column (Force character)
final_df[[gene_col]] <- as.character(final_df[[gene_col]])

plot_data <- final_df %>%
    filter(!!sym(gene_col) %in% target_genes) %>%
    select(Gene = !!sym(gene_col), LogFC = !!sym(logfc_col), PValue = !!sym(pval_col)) %>%
    mutate(Significance = case_when(
        PValue < 0.001 ~ "***",
        PValue < 0.01 ~ "**",
        PValue < 0.05 ~ "*",
        TRUE ~ "ns"
    )) %>%
    mutate(Direction = ifelse(LogFC > 0, "Upregulated (Healed)", "Downregulated (Non-healed)"))

print(plot_data)

# 4. Generate Professional Barplot
# Theme: Scientific
p <- ggplot(plot_data, aes(x = reorder(Gene, -LogFC), y = LogFC, fill = LogFC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    scale_fill_manual(
        values = c("#377EB8", "#E41A1C"),
        labels = c("Down in Healed", "Up in Healed")
    ) +
    geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
    labs(
        title = "External Validation of KIF13A Transport Module",
        subtitle = "GSE134431 (RNA-seq): Healed vs. Non-healed",
        x = "VGK Driver & Partners",
        y = "Log2 Fold Change"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none" # Simplify
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black")

# Assuming positive LogFC means Healed > Non-healed (based on KIF13A=1.65 result)
# We add an annotation arrow
p <- p + annotate("segment",
    x = 0.5, xend = 0.5, y = 0.2, yend = 1.5,
    arrow = arrow(length = unit(0.2, "cm")), color = "gray50"
) +
    annotate("text", x = 0.7, y = 1.0, label = "Healed High", angle = 90, color = "gray50")

ggsave(file.path(results_dir, "Figures/15_ExtVal_MultiGene_Validation.pdf"), p, width = 6, height = 5)

message("Advanced Figure 7A Created: 15_ExtVal_MultiGene_Validation.pdf")
