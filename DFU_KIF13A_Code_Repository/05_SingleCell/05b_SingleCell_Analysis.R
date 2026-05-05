# ==============================================================================
# Script 42: Single-Cell Mechanism Validation (Real Data Analysis)
# Purpose: Phase 4 - Validate Transport Module in Migratory Keratinocytes
# Data Source: Theocharidis et al. (2022) - GSE165816 (Real Verified Data)
# ==============================================================================

# 1. Setup & Licensing ----------------------------------------------------
rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        install.packages(pkg)
    }
}

# Core SC packages
ensure_package("Seurat")
ensure_package("tidyverse")
ensure_package("patchwork")
ensure_package("ggplot2")
ensure_package("viridis")
ensure_package("Matrix")

library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)
library(Matrix)

# 2. Data Discovery -------------------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
sc_dir <- file.path(analysis_dir, "Validation_Data/SingleCell")

message("--- Phase 4: Single Cell Analysis ---")
message("Searching for real data in: ", sc_dir)

# Pattern matching for GSE165816 raw files
files <- list.files(sc_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

if (length(files) == 0) {
    # Fallback for other formats (e.g., mtx) if tar extracted differently
    files <- list.files(sc_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
}

if (length(files) == 0) {
    stop(
        ">>> CRITICAL ERROR: No Single-Cell Data Found! \n",
        "Please run Script 91 to download GSE165816, or manually place files in Validation_Data/SingleCell.\n",
        "Analysis cannot proceed without real data."
    )
}

message(paste("Found", length(files), "data files."))

# 3. Data Loading & Merging -----------------------------------------------
# To ensure robustness and memory safety, we load a representative subset (first 5 samples)
# This is sufficient to validate the 'Cell Type' mechanism.

max_samples <- 5
files_to_load <- files[1:min(length(files), max_samples)]
message(paste("Processing first", length(files_to_load), "samples for validation..."))

obj_list <- list()

for (f in files_to_load) {
    sample_name <- gsub("_.*", "", basename(f)) # Extract GSM ID
    message(paste("Loading Sample:", sample_name))

    tryCatch(
        {
            # Handle CSV.gz
            if (grepl("csv\\.gz$", f)) {
                counts <- read.csv(gzfile(f), row.names = 1, check.names = FALSE)
                # Create Seurat Object
                obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
            }
            # Handle MTX (Future proofing)
            else {
                # Logic for MTX would go here, usually requires barcodes/features files
            }
        },
        error = function(e) {
            message("Error loading ", f, ": ", e$message)
        }
    )
}

if (length(obj_list) == 0) stop("Failed to load any Seurat objects.")

# Merge into one object
if (length(obj_list) > 1) {
    sc_obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
} else {
    sc_obj <- obj_list[[1]]
}

message("Merged Object: ", nrow(sc_obj), " genes x ", ncol(sc_obj), " cells")

# 4. Standard Pre-processing ----------------------------------------------
message("--- Pre-processing ---")
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
sc_obj <- RunUMAP(sc_obj, dims = 1:15) # Use valid dims
sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)

# 5. Mechanism Validation -------------------------------------------------
# Define Modules
transport_features <- list(c("KIF13A", "EPN1", "CLIP1", "RAB11A"))
migratory_features <- list(c("MMP9", "S100A8", "S100A9", "KRT6A"))

# Score Modules
# Using ctrl=5 to ensure it works even if random genes are sparse,
# but with real data, defaults usually work.
tryCatch(
    {
        sc_obj <- AddModuleScore(sc_obj, features = transport_features, name = "Transport_Score")
        sc_obj <- AddModuleScore(sc_obj, features = migratory_features, name = "Migration_Score")
    },
    error = function(e) {
        # Fallback with simpler control selection if binning fails
        sc_obj <- AddModuleScore(sc_obj, features = transport_features, name = "Transport_Score", ctrl = 5)
        sc_obj <- AddModuleScore(sc_obj, features = migratory_features, name = "Migration_Score", ctrl = 5)
    }
)

# Rename
sc_obj$Transport_Module <- sc_obj$Transport_Score1
sc_obj$Migration_State <- sc_obj$Migration_Score1

# 6. Visualization --------------------------------------------------------
out_pdf <- file.path(analysis_dir, "Figures/18_SingleCell_Validation.pdf")
pdf(out_pdf, width = 12, height = 8)

# Plots
p1 <- FeaturePlot(sc_obj, features = "Transport_Module", pt.size = 0.5, order = TRUE) +
    scale_color_viridis(option = "magma") + ggtitle("Transport Network Activation")

p2 <- FeaturePlot(sc_obj, features = "Migration_State", pt.size = 0.5, order = TRUE) +
    scale_color_viridis(option = "cividis") + ggtitle("Migratory Phenotype")

p3 <- DimPlot(sc_obj, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")

# Correlation Plot (Transport vs Migration)
df_cor <- data.frame(
    Transport = sc_obj$Transport_Module,
    Migration = sc_obj$Migration_State,
    Cluster = sc_obj$seurat_clusters
)
cor_val <- cor(df_cor$Transport, df_cor$Migration, method = "pearson")

p4 <- ggplot(df_cor, aes(x = Migration, y = Transport)) +
    geom_point(alpha = 0.1, color = "grey") +
    geom_smooth(method = "lm", color = "red") +
    theme_bw() +
    ggtitle(paste0("Correlation: R = ", round(cor_val, 2))) +
    labs(x = "Migratory State Score", y = "Transport Module Score")

# Layout
layout_design <- "
AABB
CCDD
"
print(p1 + p2 + p3 + p4 + plot_layout(design = layout_design))

dev.off()
message("Analysis Complete. Figure saved to: ", out_pdf)
