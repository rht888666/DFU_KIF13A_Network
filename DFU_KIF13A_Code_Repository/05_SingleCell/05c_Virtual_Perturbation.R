# ==============================================================================
# Script 44: Single-Cell Virtual Perturbation (Virtual OE & KO)
# Purpose: Deep Mechanism Digging using "In Silico Perturbation"
# Techniques:
#   1. Virtual Knockout (Network Topology Collapse Simulation)
#   2. Virtual Overexpression (Digital Twin / High-Low Comparison)
#   3. Rescue Simulation (Can KIF13A restoration fix Non-healed cells?)
# Data Source: GSE165816 (Real Verified Data)
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

# Core Packages
ensure_package("Seurat")
ensure_package("tidyverse")
ensure_package("igraph")
ensure_package("patchwork")
ensure_package("viridis")
ensure_package("Matrix")
ensure_package("ggrepel")

# Robustness Options
options(lifecycle_verbosity = "warning")
options(future.globals.maxSize = 8000 * 1024^2)

library(Seurat)
library(tidyverse)
library(igraph)
library(patchwork)
library(viridis)
library(ggrepel)

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
library(Matrix)

# 2. Data Loading (Fast Track) --------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
sc_dir <- file.path(analysis_dir, "Validation_Data/SingleCell")

# Try to find processed object first (if saved), else reload raw
# For this script, we assume we need to reload as per Script 42 logic
files <- list.files(sc_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) files <- list.files(sc_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)

if (length(files) == 0) stop("No Single-Cell Data Found! Please check Script 91/42.")

# Load a representative subset for analysis (First 5 samples to save memory, or all if feasible)
# We need enough cells to build a robust network.
max_samples <- 8 # Increased to get enough Non-healed cells
files_to_load <- files[1:min(length(files), max_samples)]

obj_list <- list()
message("Loading real single-cell data for Virtual Perturbation...")
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

# Standard Preprocessing
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & percent.mt < 15)
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, nfeatures = 2000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, verbose = FALSE)
sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)
sc_obj <- RunUMAP(sc_obj, dims = 1:15)

message("Data Ready: ", nrow(sc_obj), " genes x ", ncol(sc_obj), " cells")

# 3. Identify Target Cell Population (Keratinocytes) -----------------------
# We look for clusters high in KRT14, KRT10, or KRT6A
# Robust Cluster Identification
markers <- c("KRT14", "KRT10", "KRT6A")
valid_markers <- markers[markers %in% rownames(sc_obj)]

if (length(valid_markers) == 0) {
    message("Warning: Key Keratinocyte markers not found. Using Cluster 0 as default.")
    target_cluster <- 0
} else {
    sc_obj <- AddModuleScore(sc_obj, features = list(valid_markers), name = "Keratinocyte_Score", nbin = 10)
    # Check if column added
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

# 4. PART A: Single-Cell Virtual Knockout (The 'Collapse' Simulation) -----
message("--- Executing Virtual Knockout of KIF13A ---")

# We implement a graph-theoretic perturbation
# 1. Build Co-expression Network (Adjacency Matrix) on Top Variable Genes
top_genes <- VariableFeatures(kc_obj)[1:1000] # Top 1000 for speed/robustness
# Ensure KIF13A and key Cargo are in list
targets <- c("KIF13A", "MMP9", "ITGB1", "VEGFA", "S100A8")
top_genes <- unique(c(top_genes, targets))
top_genes <- top_genes[top_genes %in% rownames(kc_obj)]

# Step A: Connectivity Map (Co-expression)
# Using `get_data` for robust access
expr_mat <- get_data(kc_obj)[top_genes, ]
expr_mat <- as.matrix(expr_mat)


# Calculate Correlation Matrix (WGCNA style)
cor_mat <- cor(t(expr_mat), method = "pearson")
cor_mat[is.na(cor_mat)] <- 0

# Soft Thresholding to create Adjacency
beta <- 6 # Typical power
adj_mat <- abs(cor_mat)^beta
diag(adj_mat) <- 0

# Baseline Centrality (Eigenvector)
# We use PageRank as it's robust for directed flow simulations, or Eigenvector
g_base <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
pr_base <- page_rank(g_base)$vector

# The Perturbation: Remove KIF13A Node
# Function to simulate KO
simulate_ko <- function(adj, gene) {
    if (!gene %in% rownames(adj)) {
        return(NULL)
    }
    adj_ko <- adj
    adj_ko[gene, ] <- 0
    adj_ko[, gene] <- 0
    g_ko <- graph_from_adjacency_matrix(adj_ko, mode = "undirected", weighted = TRUE)
    pr_ko <- page_rank(g_ko)$vector
    return(pr_ko)
}

# Run KO
if ("KIF13A" %in% rownames(adj_mat)) {
    pr_ko <- simulate_ko(adj_mat, "KIF13A")

    # Calculate Impact: % Change in Centrality
    # "Who loses importance when KIF13A is gone?"
    impact <- (pr_base - pr_ko) / pr_base
    impact_df <- data.frame(
        Gene = names(impact),
        Centrality_Base = pr_base,
        Centrality_KO = pr_ko,
        Impact_Score = impact
    )

    # Remove KIF13A itself from result
    impact_df <- impact_df[impact_df$Gene != "KIF13A", ]
    impact_df <- impact_df %>% arrange(desc(Impact_Score))

    message("Top genes impacted by KIF13A Virtual KO:")
    print(head(impact_df, 10))

    # Save result
    write.csv(impact_df, file.path(analysis_dir, "Tables/44_SC_Virtual_KO_Results.csv"))
} else {
    message("KIF13A not sufficiently expressed in this subset to build network.")
    impact_df <- NULL
}

# 5. PART B: Single-Cell Virtual Overexpression (The 'Rescue' Simulation) --
message("--- Executing Virtual Overexpression (In Silico Perturbation) ---")

# Strategy: "Digital Twin" / "Natural Perturbation"
# We define "Overexpression" not by adding fake counts, but by identifying
# cells that SPONTANEOUSLY overexpress KIF13A within the Non-Healed-like population
# and asking: "Do they perform better?"

# 1. Define 'State'
# We need to distinguish Healed vs Non-healed cells.
# Since metadata might be tricky in merged objects, we use gene signatures.
# Non-Healing Signature: S100A8, S100A9, MMP9 High
# Healing Signature: MKI67 (proliferation), COL1A1 (remodeling) - simplified

# Let's assume the dataset contains mix of conditions.
# We create a 'Virtual OE' experiment:
# Group 1: KIF13A High Expression (Top 20%)
# Group 2: KIF13A Low Expression (Bottom 20%)
# Control: Cells must be from the same cluster (Keratinocytes) to avoid confounding.

kif_vals <- tryCatch(
    {
        get_data(kc_obj)["KIF13A", ]
    },
    error = function(e) {
        # If KIF13A is missing or error
        rep(0, ncol(kc_obj))
    }
)
high_cutoff <- quantile(kif_vals[kif_vals > 0], 0.75) # Top 25% of expressors
low_cutoff <- 0 # Strictly zero or very low

kc_obj$KIF13A_Status <- case_when(
    kif_vals >= high_cutoff ~ "Virtual_OE",
    kif_vals == 0 ~ "Virtual_KO",
    TRUE ~ "Neutral"
)

# 2. Analyze the 'Rescue' Effect
# Does Virtual_OE status drive 'Migration' or 'Differentiation'?
# Define functional scores
mig_genes <- list(c("MMP1", "MMP9", "PLAU", "ITGB1")) # Migration
inflam_genes <- list(c("S100A8", "S100A9", "IL1B", "CXCL8")) # Inflammation

kc_obj <- AddModuleScore(kc_obj, features = mig_genes, name = "Migration_Score", nbin = 10)
kc_obj <- AddModuleScore(kc_obj, features = inflam_genes, name = "Inflammation_Score", nbin = 10)

# Rename module columns
kc_obj$Migration <- kc_obj$Migration_Score1
kc_obj$Inflammation <- kc_obj$Inflammation_Score1

# Statistical Test (Likelihood of Improved Phenotype)
# Compare 'Migration' in OE vs KO groups
test_res <- wilcox.test(
    kc_obj$Migration[kc_obj$KIF13A_Status == "Virtual_OE"],
    kc_obj$Migration[kc_obj$KIF13A_Status == "Virtual_KO"]
)

message("Virtual Overexpression Test (Migration): P-value = ", test_res$p.value)

# 6. Visualization --------------------------------------------------------
pdf(file.path(analysis_dir, "Figures/44_Virtual_Perturbation_Analysis_v2.pdf"), width = 12, height = 10)

layout_design <- "
AABB
CCDD
"

# Plot A: Virtual KO Impact (Volcano-like)
if (!is.null(impact_df)) {
    # Highlight top genes
    top_hits <- head(impact_df$Gene, 10)
    p1 <- ggplot(impact_df, aes(x = Centrality_Base, y = Impact_Score)) +
        geom_point(color = "grey", alpha = 0.6) +
        geom_point(data = subset(impact_df, Gene %in% top_hits), color = "red", size = 2) +
        geom_text_repel(data = subset(impact_df, Gene %in% top_hits), aes(label = Gene), size = 3, max.overlaps = 20) +
        theme_minimal() +
        labs(
            title = "Virtual Knockout Impact (Network Collapse)",
            x = "Baseline Centrality", y = "Centrality Loss after KIF13A Removal",
            subtitle = "Top genes = 'Stranded Cargo'"
        )
} else {
    p1 <- ggplot() +
        annotate("text", x = 1, y = 1, label = "Not enough KIF13A data")
}

# Plot B: KIF13A Expression Distribution (Defining the groups)
p2 <- VlnPlot(kc_obj, features = "KIF13A", group.by = "seurat_clusters", pt.size = 0) +
    ggtitle("KIF13A Expression in Keratinocytes") + theme(legend.position = "none")

# Plot C: Virtual Overexpression Effect (The Rescue)
# Boxplot of Migration Score by Virtual Status
p3 <- ggplot(
    kc_obj@meta.data %>% filter(KIF13A_Status != "Neutral"),
    aes(x = KIF13A_Status, y = Migration, fill = KIF13A_Status)
) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.2, fill = "white") +
    scale_fill_manual(values = c("Virtual_KO" = "#E69F00", "Virtual_OE" = "#56B4E9")) +
    theme_bw() +
    labs(
        title = "Virtual Overexpression Simulation",
        subtitle = paste0("Does KIF13A OE drive migration? P = ", format.pval(test_res$p.value, digits = 3)),
        y = "Migratory Phenotype Score", x = "Simulated Status"
    )

# Plot D: Correlation Network of Top Regulators (Virtual OE Context)
# Show what KIF13A correlates with *specifically* in the High group
# (Simplified to Scatter plot for robustness)
p4 <- FeatureScatter(kc_obj, feature1 = "KIF13A", feature2 = "Migration", group.by = "KIF13A_Status") +
    ggtitle("KIF13A vs Migration Correlation") +
    scale_color_manual(values = c("Virtual_KO" = "lightgrey", "Neutral" = "grey", "Virtual_OE" = "blue"))

print(p1 + p2 + p3 + p4 + plot_layout(design = layout_design))

dev.off()

message("Virtual Perturbation Analysis Complete.")
message("Results saved to: Figures/44_Virtual_Perturbation_Analysis_v2.pdf")
