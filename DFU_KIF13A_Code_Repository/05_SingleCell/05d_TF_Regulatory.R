# ==============================================================================
# Script 45: Correlation-Based TF Inference (Simplified Regulatory Network)
# Purpose: Identify Potential Upstream Regulators of KIF13A
# Method: Calculate Spearman correlation between KIF13A and all known Human TFs
# Goal: Find NEGATIVE regulators (Repressors) linked to inflammation
# Note: This is a lightweight alternative to full SCENIC analysis.
# ==============================================================================

# 1. Setup ----------------------------------------------------------------
rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Robustness Options
options(lifecycle_verbosity = "warning")
options(future.globals.maxSize = 8000 * 1024^2)

library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)

# Helper function for Seurat v5 compatibility
get_data <- function(obj) {
    # Seurat v5/v4 compatibility wrapper
    if (packageVersion("Seurat") >= "5.0.0") {
        # Auto-join layers if 'data' layer is missing (common in v5 split objects)
        if (!"data" %in% Layers(obj)) {
            obj <- JoinLayers(obj)
        }
        return(GetAssayData(obj, layer = "data"))
    } else {
        return(GetAssayData(obj, slot = "data"))
    }
}
# 2. Data Loading (Fast Track) --------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
sc_dir <- file.path(analysis_dir, "Validation_Data/SingleCell")

message("Loading Single-Cell Data for TF Correlation...")

files <- list.files(sc_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) files <- list.files(sc_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No Single-Cell Data Found! Please check Script 91/42.")

max_samples <- 5
files_to_load <- files[1:min(length(files), max_samples)]

obj_list <- list()
for (f in files_to_load) {
    sample_name <- gsub("_.*", "", basename(f))
    tryCatch(
        {
            if (grepl("csv\\.gz$", f)) {
                counts <- read.csv(gzfile(f), row.names = 1, check.names = FALSE)
                obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
            }
        },
        error = function(e) message("Skip ", f)
    )
}

if (length(obj_list) == 0) stop("Data Load Failed.")
sc_obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))

# Basic Preprocessing
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & percent.mt < 15)
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, nfeatures = 2000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, verbose = FALSE)
sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)
sc_obj <- RunUMAP(sc_obj, dims = 1:15)

message("Data loaded: ", ncol(sc_obj), " cells.")

# 3. Focus on Keratinocytes ----------------------------------------------
# Robust Cluster Identification
markers <- c("KRT14", "KRT10")
valid_markers <- markers[markers %in% rownames(sc_obj)]

if (length(valid_markers) == 0) {
    message("Warning: Key Keratinocyte markers not found. Using Cluster 0 as default.")
    target_cluster <- 0
} else {
    sc_obj <- AddModuleScore(sc_obj, features = list(valid_markers), name = "Keratinocyte_Score", nbin = 10)
    if ("Keratinocyte_Score1" %in% colnames(sc_obj@meta.data)) {
        cluster_scores <- aggregate(sc_obj$Keratinocyte_Score1, by = list(sc_obj$seurat_clusters), median)
        target_cluster <- cluster_scores$Group.1[which.max(cluster_scores$x)]
    } else {
        message("AddModuleScore failed. Using Cluster 0.")
        target_cluster <- 0
    }
}

kc_obj <- subset(sc_obj, idents = target_cluster)
# Seurat v5: Ensure layers are joined for downstream matrix operations
if (packageVersion("Seurat") >= "5.0.0") {
    kc_obj <- JoinLayers(kc_obj)
}
message("Subsetting to ", ncol(kc_obj), " Keratinocytes.")

# 4. Correlation Analysis ------------------------------------------------
message("--- Analyzing TF Correlations ---")

# Step A: Get Human TF List
# Since we don't have the database installed, we define a core list of known TFs manually
# This covers major families: NFkB, AP-1, STAT, KLF, FOX, SOX, etc.
core_tfs <- c(
    "RELA", "NFKB1", "NFKB2", "REL", "RELB", # NFkB Family (Inflammation)
    "STAT1", "STAT3", "STAT5A", "STAT6", # STAT Family (Inflammation)
    "JUN", "JUNB", "JUND", "FOS", "FOSB", # AP-1 Family (Stress/Migration)
    "KLF4", "KLF5", "KLF2", # KLF Family (Differentiation)
    "TP63", "TP53", # p53/p63 (Cell Cycle/Stemness)
    "MYC", "HIF1A", "EP300", "CREB1", # General Regulators
    "SNAI1", "SNAI2", "TWIST1", "ZEB1", # EMT (Migration)
    "GRHL1", "GRHL2", "GRHL3", # Grainyhead (Epithelial)
    "OVOL1", "OVOL2", # Ovol (Epithelial)
    "IRF1", "IRF3", "IRF7", # Interferon (Inflammation)
    "ETV4", "ETV5", # PEA3 (Migration)
    "ELF3", "EHF" # ESE (Epithelial)
)

# Filter for TFs present in dataset
available_tfs <- intersect(rownames(kc_obj), core_tfs)
message("Analyzing ", length(available_tfs), " Transcription Factors.")

if (!"KIF13A" %in% rownames(kc_obj)) stop("KIF13A not expressed in subset.")

# Step B: Calculate Correlations
# Get Expression Matrix (Genes x Cells)
# Use 'data' slot# Step B: Calculate Correlations
# Robust data access
expr_mat <- get_data(kc_obj)[c("KIF13A", available_tfs), ]
expr_mat <- as.matrix(expr_mat)

# Spearman Correlation (KIF13A vs TFs)
# Spearman Correlation (KIF13A vs TFs)
kif13a_expr <- expr_mat["KIF13A", ]
tf_expr <- expr_mat[available_tfs, , drop = FALSE]

cor_res <- apply(tf_expr, 1, function(x) {
    cor(x, kif13a_expr, method = "spearman")
})

cor_df <- data.frame(TF = names(cor_res), Correlation = cor_res)
cor_df <- cor_df %>% arrange(desc(Correlation))

# Step C: Identify Candidates
message("Top POSITIVE Correlators (Potential Activators):")
print(head(cor_df, 8))

message("Top NEGATIVE Correlators (Potential Repressors):")
print(tail(cor_df, 8))

# Save
write.csv(cor_df, file.path(analysis_dir, "Tables/45_TF_Correlation_Results.csv"))

# 5. Visualization -------------------------------------------------------
pdf(file.path(analysis_dir, "Figures/45_Regulatory_QuickLook.pdf"), width = 10, height = 6)

# Plot A: Barplot
cor_df$Type <- ifelse(cor_df$Correlation > 0, "Pos", "Neg")
cor_df$TF <- factor(cor_df$TF, levels = cor_df$TF[order(cor_df$Correlation)])

p1 <- ggplot(cor_df, aes(x = TF, y = Correlation, fill = Type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Pos" = "forestgreen", "Neg" = "firebrick")) +
    theme_bw() +
    labs(
        title = "KIF13A Regulatory Landscape",
        subtitle = "Correlation with Key TFs in Keratinocytes",
        y = "Spearman Correlation", x = "Transcription Factor"
    )

# Plot B: Scatter of Top Negative TF
top_neg_tf <- tail(cor_df$TF, 1)
p2 <- FeatureScatter(kc_obj, feature1 = as.character(top_neg_tf), feature2 = "KIF13A") +
    ggtitle(paste0("Potential Repressor: ", top_neg_tf)) +
    theme(legend.position = "none")

print(p1 + p2)
dev.off()

message("Quick Regulatory Analysis Complete.")
message("Results saved to: Figures/45_Regulatory_QuickLook.pdf")
