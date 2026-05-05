# ============================================================================
# DFU Real Data Analysis - Single Cell Validation
# Script 05: KIF13A Single-Cell Localization (GSE165816)
# Target: Nature Medicine / JAMA / Lancet
# Note: STRICTLY REAL DATA. Downloads data if not present.
# ============================================================================

rm(list = ls())
gc()

# ============================================================================
# 1. 环境设置
# ============================================================================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 3600) # 增加下载超时时间

suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
data_dir <- file.path(results_dir, "Data/SingleCell")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "KIF13A"

# ============================================================================
# 2. 检�?获取单细胞数�?(GSE165816)
# ============================================================================

message("Checking for Single-Cell Data (GSE165816)...")

# 由于完整数据集非常大，为了演示真实性，我们尝试下载处理好的表达矩阵或其子集
# 这里我们假设用户可能没有下载数GB的原始数据�?
# 为了代码的可执行性，我们将检查本地是否有数据，如果没有，
# **我们将尝试使�?SeuratData (如果适用) 或者提示必须手动下�?*
# 但为了满�?全自�?的要求，我们将尝试下载一个轻量级的替代方案（如果存在），
# 或者必须依赖用户提供文件�?

# 鉴于环境限制，我们先检查本地特定路�?
local_file <- file.path(data_dir, "GSE165816_processed.rds")

if (file.exists(local_file)) {
    message("Loading local processed Seurat object...")
    sc_obj <- readRDS(local_file)
} else {
    message("Local processed file not found.")
    message("Attempting to download GSE165816 metadata and matrix...")

    # 真实情况：我们需要从 GEO 下载�?
    # 这里我们模拟下载逻辑：如果检测不到文件，我们将停止并报错�?
    # 因为在无交互环境下下�?3GB+ 数据是不现实的�?
    # *除非* 我们能找到一个精简版�?

    # 为了保证流程走通，我将检查是否可以从项目旧目录找到数�?
    old_data_path <- "E:/人工智能学习/任海涛虚拟基因敲�?Data/GSE165816"
    if (dir.exists(old_data_path)) {
        message("Found data in old project directory. Importing...")
        # 假设�?10X 格式
        tryCatch(
            {
                sc_data <- Read10X(data.dir = old_data_path)
                sc_obj <- CreateSeuratObject(counts = sc_data, project = "DFU_Real")
            },
            error = function(e) {
                message("Error reading 10X data: ", e$message)
                stop("Please ensure GSE165816 10X files (barcodes, features, matrix) are in 'Data/SingleCell'")
            }
        )
    } else {
        # 既然我们承诺�?真实数据"，如果真的没有文件，
        # 我们不能生成"模拟"的�?
        # 我们必须生成一�?�?报告，告知用户缺数据�?

        sink(file.path(results_dir, "Reports/05_SingleCell_Status.txt"))
        cat("================================================\n")
        cat("          SINGLE CELL DATA MISSING              \n")
        cat("================================================\n")
        cat("CRITICAL: Real single-cell data (GSE165816) not found.\n")
        cat("Action Required: Please download GSE165816_RAW.tar from GEO\n")
        cat("and extract to:", data_dir, "\n")
        sink()

        stop("REAL DATA MISSING: Cannot perform Single-Cell analysis without GSE165816 file. Stopping to avoid simulation.")
    }
}

# ============================================================================
# 3. 数据处理 (如果尚未处理)
# ============================================================================

if (!"pca" %in% names(sc_obj@reductions)) {
    message("Preprocessing Seurat object...")
    sc_obj <- NormalizeData(sc_obj)
    sc_obj <- FindVariableFeatures(sc_obj)
    sc_obj <- ScaleData(sc_obj)
    sc_obj <- RunPCA(sc_obj)
    sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
    sc_obj <- FindClusters(sc_obj, resolution = 0.5)
    sc_obj <- RunUMAP(sc_obj, dims = 1:15)
}

# ============================================================================
# 4. KIF13A 定位分析
# ============================================================================

message(paste("Analyzing", target_gene, "expression..."))

# 检查基因是否存�?
if (!target_gene %in% rownames(sc_obj)) {
    stop(paste(target_gene, "not found in single-cell matrix!"))
}

# 1. FeaturePlot (UMAP)
p1 <- FeaturePlot(sc_obj, features = target_gene, pt.size = 0.5) +
    ggtitle(paste(target_gene, "Expression in DFU Tissue"))

# 2. VlnPlot (按簇/细胞类型)
# 如果有细胞类型注释，使用之；否则使用 Cluster
group_col <- if ("cell_type" %in% colnames(sc_obj@meta.data)) "cell_type" else "seurat_clusters"
p2 <- VlnPlot(sc_obj, features = target_gene, group.by = group_col) + NoLegend()

# 3. DotPlot (邻居基因共定�?
neighbors <- c(target_gene, "EPN1", "ANXA1", "ZNF185") # 来自 Script 04 的发�?
valid_neighbors <- neighbors[neighbors %in% rownames(sc_obj)]
p3 <- DotPlot(sc_obj, features = valid_neighbors, group.by = group_col) +
    RotatedAxis() +
    ggtitle("Co-expression of KIF13A and Neighbors")

# 组合�?
p_combined <- (p1 | p2) / p3
ggsave(file.path(results_dir, "Figures/05_KIF13A_SingleCell_Real.pdf"), p_combined, width = 12, height = 10)

# ============================================================================
# 5. 生成报告
# ============================================================================

sink(file.path(results_dir, "Reports/05_SingleCell_Findings.txt"))
cat("================================================\n")
cat("          SINGLE CELL VALIDATION (REAL)         \n")
cat("================================================\n")
cat("Target Gene:", target_gene, "\n")
cat("Data Source: GSE165816 (Real)\n\n")

# 计算表达统计
expr_vals <- GetAssayData(sc_obj, slot = "data")[target_gene, ]
pos_cells <- sum(expr_vals > 0)
total_cells <- length(expr_vals)
cat("Expression Stats:\n")
cat("   - Positive Cells:", pos_cells, "(", round(pos_cells / total_cells * 100, 2), "%)\n")
cat("   - Max Expression:", max(expr_vals), "\n\n")

# 按簇统计
cluster_stats <- data.frame(
    Cluster = sc_obj@meta.data[[group_col]],
    Expr = expr_vals
) %>%
    group_by(Cluster) %>%
    summarise(
        Mean_Expr = mean(Expr),
        Pct_Pos = sum(Expr > 0) / n() * 100
    ) %>%
    arrange(desc(Mean_Expr))

cat("Top Expressing Cell Types/Clusters:\n")
print(head(cluster_stats, 5))

sink()

message("Script 05 Completed. Real Single-Cell validation done.")
