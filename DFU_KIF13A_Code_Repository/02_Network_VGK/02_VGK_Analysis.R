# ============================================================================
# DFU Virtual Gene Knockout - Real Data Analysis Suite (Authentic Mode)
# Script 02: Real VGK Analysis - Identifying Key Drivers
# Target: Nature Medicine / JAMA / Lancet
# Note: STRICTLY REAL DATA ONLY.
# ============================================================================

# 清空工作环境
rm(list = ls())
gc()

# ============================================================================
# 1. 环境设置
# ============================================================================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

suppressPackageStartupMessages({
    library(tidyverse)
    library(igraph)
    library(foreach)
    library(doParallel)
})

# 设置并行计算
n_cores <- min(parallel::detectCores() - 1, 4)
registerDoParallel(n_cores)

# ============================================================================
# 2. 加载真实数据
# ============================================================================

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
data_file <- file.path(results_dir, "Data/01_REAL_processed_data.RData")

if (!file.exists(data_file)) stop("CRITICAL ERROR: Real processed data not found! Run Script 01 first.")

load(data_file)
message("Loaded REAL processed data.")

# ============================================================================
# 3. 构建网络函数 (VGK Core)
# ============================================================================

# 网络构建函数
build_network <- function(expr_matrix, threshold = 0.5) {
    # 计算相关�?(Pearson)
    cor_mat <- cor(t(expr_matrix), method = "pearson", use = "pairwise.complete.obs")
    # 构建邻接矩阵
    adj_mat <- abs(cor_mat)
    adj_mat[adj_mat < threshold] <- 0
    diag(adj_mat) <- 0
    # 创建igraph对象
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
    return(g)
}

# 网络凝聚力计算函�?(Cohesion Metrics)
calc_cohesion <- function(g) {
    if (vcount(g) == 0 || ecount(g) == 0) {
        return(0)
    }

    # 指标1: 密度 (Density)
    density <- edge_density(g)

    # 指标2: 全局聚类系数 (Transitivity)
    clustering <- transitivity(g, type = "global")
    if (is.nan(clustering)) clustering <- 0

    # 指标3: 平均路径长度的倒数 (Efficiency) - 计算量大，暂用最大连通分量占比替�?
    # efficiency <- 1 / mean_distance(g)

    # 综合凝聚力指�?(CI) = (Density + Clustering) / 2
    return((density + clustering) / 2)
}

# VGK 扰动函数
perform_vgk <- function(g, gene_list) {
    original_score <- calc_cohesion(g)

    impact_scores <- foreach(gene = gene_list, .combine = rbind, .packages = "igraph", .export = "calc_cohesion") %dopar% {
        if (gene %in% V(g)$name) {
            # 虚拟敲除: 移除节点及其连边
            g_knockout <- delete_vertices(g, gene)
            new_score <- calc_cohesion(g_knockout)
            # 影响�?= 原始凝聚�?- 敲除后凝聚力
            impact <- original_score - new_score
            return(data.frame(gene = gene, impact = impact))
        } else {
            return(data.frame(gene = gene, impact = 0))
        }
    }
    return(impact_scores)
}

# ============================================================================
# 4. 执行真实 VGK 分析
# ============================================================================

message("Starting Real VGK Analysis...")

# 筛选高变基�?(Top 2000) 以构建核心调控网�?
# 全基因组网络计算量太大，且包含大量噪�?
gene_vars <- apply(expr_scaled, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]

# 确保之前的候选基�?(�?ZSWIM8) 在列表中，以便验�?
candidates <- c("ZSWIM8", "RAB3IL1", "GNB2", "AKR1B1", "MMP1", "COL1A1")
top_genes <- unique(c(top_genes, candidates))
top_genes <- top_genes[top_genes %in% rownames(expr_scaled)]

message(paste("Analyzing network with", length(top_genes), "genes."))

# 分组样本
healed_samples <- metadata$sample[metadata$healed == "Yes"]
nonhealed_samples <- metadata$sample[metadata$healed == "No"]

message(paste("Healed samples:", length(healed_samples)))
message(paste("Non-healed samples:", length(nonhealed_samples)))

# 构建愈合网络 (Healed Network)
message("Building Healed Network...")
expr_healed <- expr_scaled[top_genes, healed_samples]
g_healed <- build_network(expr_healed, threshold = 0.6) # 使用较高阈值确保网络稳健�?

# 构建不愈合网�?(Non-healed Network)
message("Building Non-healed Network...")
expr_nonhealed <- expr_scaled[top_genes, nonhealed_samples]
g_nonhealed <- build_network(expr_nonhealed, threshold = 0.6)

# 执行 VGK
message("Performing Virtual Knockout on Healed Network...")
res_healed <- perform_vgk(g_healed, top_genes)
colnames(res_healed)[2] <- "Impact_Healed"

message("Performing Virtual Knockout on Non-healed Network...")
res_nonhealed <- perform_vgk(g_nonhealed, top_genes)
colnames(res_nonhealed)[2] <- "Impact_Nonhealed"

# ============================================================================
# 5. 结果整合与差异影响力分析 (DI Score)
# ============================================================================

message("Calculating Differential Impact (DI) Scores...")

vgk_results <- merge(res_healed, res_nonhealed, by = "gene")

# DI Score = Impact_Nonhealed - Impact_Healed
# 正�? 维持不愈合网络所必需 (不愈合特异性驱动因�?
# 负�? 维持愈合网络所必需 (愈合特异性驱动因�?
vgk_results$DI_Score <- vgk_results$Impact_Nonhealed - vgk_results$Impact_Healed

# 排序
vgk_results <- vgk_results %>% arrange(desc(DI_Score))

# ============================================================================
# 6. 保存结果
# ============================================================================

message("Saving VGK Results...")

write.csv(vgk_results, file.path(results_dir, "Tables/02_Real_VGK_Results.csv"), row.names = FALSE)

# 保存 Top 50 基因列表
top_drivers <- head(vgk_results, 50)
write.csv(top_drivers, file.path(results_dir, "Tables/02_Real_Top_Drivers.csv"), row.names = FALSE)

# 生成简�?
sink(file.path(results_dir, "Reports/02_VGK_Summary.txt"))
cat("================================================\n")
cat("          REAL VGK ANALYSIS SUMMARY             \n")
cat("================================================\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("1. Network Parameters:\n")
cat("   - Gene Universe:", length(top_genes), "\n")
cat("   - Correlation Threshold: 0.6\n\n")
cat("2. Top 10 Non-Healing Drivers (Positive DI):\n")
print(head(vgk_results[, c("gene", "DI_Score")], 10))
cat("\n3. Top 10 Healing Drivers (Negative DI):\n")
print(tail(vgk_results[, c("gene", "DI_Score")], 10))
cat("\n4. ZSWIM8 Status:\n")
zswim8_res <- vgk_results[vgk_results$gene == "ZSWIM8", ]
if (nrow(zswim8_res) > 0) {
    print(zswim8_res)
    rank <- which(vgk_results$gene == "ZSWIM8")
    cat("   Rank:", rank, "/", nrow(vgk_results), "\n")
} else {
    cat("   ZSWIM8 not found in top variable genes.\n")
}
sink()

message("Script 02 Completed. Check Reports/02_VGK_Summary.txt for findings.")
