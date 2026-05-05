# ============================================================================
# DFU Real Data Analysis - Mechanism Exploration
# Script 04: KIF13A Functional Characterization & Mechanism Decoding
# Target: Nature Medicine / JAMA / Lancet
# Note: STRICTLY REAL DATA ONLY.
# ============================================================================

rm(list = ls())
gc()

# ============================================================================
# 1. 环境设置
# ============================================================================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

suppressPackageStartupMessages({
    library(tidyverse)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(cowplot)
})

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
load(file.path(results_dir, "Data/01_REAL_processed_data.RData"))

target_gene <- "KIF13A"

# ============================================================================
# 2. 识别 KIF13A 共表达邻�?(Co-expression Neighbors)
# ============================================================================

message(paste("Identifying co-expression neighbors for", target_gene, "..."))

if (!target_gene %in% rownames(expr_log)) stop(paste("Gene", target_gene, "not found in expression data!"))

# 计算 KIF13A 与所有基因的相关�?
target_expr <- as.numeric(expr_log[target_gene, ])
cor_results <- apply(expr_log, 1, function(x) cor(x, target_expr, method = "pearson"))

# 筛选显著相关的基因 (Top 300, 且相关�?> 0.4)
cor_df <- data.frame(Gene = names(cor_results), Correlation = cor_results) %>%
    filter(Gene != target_gene) %>%
    arrange(desc(abs(Correlation)))

top_neighbors <- cor_df %>%
    filter(abs(Correlation) > 0.4) %>%
    head(300)
message(paste("Found", nrow(top_neighbors), "significant neighbors for", target_gene))

write.csv(top_neighbors, file.path(results_dir, "Tables/04_KIF13A_Neighbors.csv"), row.names = FALSE)

# ============================================================================
# 3. 功能富集分析 (GO & KEGG)
# ============================================================================

message("Performing functional enrichment analysis...")

# 转换基因ID
gene_list <- bitr(top_neighbors$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if (nrow(gene_list) > 0) {
    # GO BP Enrichment
    ego <- enrichGO(
        gene = gene_list$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    # KEGG Enrichment
    kk <- enrichKEGG(
        gene = gene_list$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.05
    )

    # 保存结果
    if (!is.null(ego)) write.csv(as.data.frame(ego), file.path(results_dir, "Tables/04_KIF13A_GO_Enrichment.csv"))
    if (!is.null(kk)) write.csv(as.data.frame(kk), file.path(results_dir, "Tables/04_KIF13A_KEGG_Enrichment.csv"))

    # 可视�?Top Pathways
    if (!is.null(ego) && nrow(ego) > 0) {
        p_go <- dotplot(ego, showCategory = 15) +
            ggtitle(paste("GO Enrichment of", target_gene, "Neighbors")) +
            theme(plot.title = element_text(size = 12, face = "bold"))
        ggsave(file.path(results_dir, "Figures/04_KIF13A_GO_Dotplot.pdf"), p_go, width = 10, height = 8)
    }

    if (!is.null(kk) && nrow(kk) > 0) {
        p_kegg <- dotplot(kk, showCategory = 15) +
            ggtitle(paste("KEGG Enrichment of", target_gene, "Neighbors"))
        ggsave(file.path(results_dir, "Figures/04_KIF13A_KEGG_Dotplot.pdf"), p_kegg, width = 10, height = 8)
    }
} else {
    message("No valid Entrez IDs found for enrichment.")
}

# ============================================================================
# 4. 关键通路相关性验�?(Angiogenesis & Vesicle Transport)
# ============================================================================

message("Verifying correlations with specific pathways...")

# 定义通路标记基因
pathway_markers <- list(
    "Angiogenesis" = c("VEGFA", "KDR", "FLT1", "PECAM1", "VWF"),
    "Vesicle_Transport" = c("RAB5A", "RAB7A", "STX1A", "VAMP2", "SNAP25"),
    "Inflammation" = c("IL6", "TNF", "CXCL8", "IL1B", "CCL2"),
    "ECM_Remodeling" = c("COL1A1", "COL3A1", "MMP9", "MMP2", "TIMP1")
)

plot_list <- list()

for (pathway in names(pathway_markers)) {
    markers <- pathway_markers[[pathway]]
    valid_markers <- markers[markers %in% rownames(expr_log)]

    if (length(valid_markers) > 0) {
        # 计算通路评分 (平均表达�?
        pathway_score <- colMeans(expr_log[valid_markers, , drop = FALSE])

        plot_data <- data.frame(
            KIF13A = target_expr,
            Pathway_Score = pathway_score,
            Group = metadata$healed
        )

        # 计算相关�?
        cor_val <- cor(plot_data$KIF13A, plot_data$Pathway_Score, method = "pearson")
        p_val <- cor.test(plot_data$KIF13A, plot_data$Pathway_Score)$p.value

        p <- ggplot(plot_data, aes(x = KIF13A, y = Pathway_Score)) +
            geom_point(aes(color = Group), alpha = 0.6) +
            geom_smooth(method = "lm", color = "black", se = TRUE) +
            scale_color_manual(values = c("Yes" = "#4DAF4A", "No" = "#E41A1C")) +
            labs(
                title = paste(target_gene, "vs", pathway),
                subtitle = paste("R =", round(cor_val, 3), ", P =", format.pval(p_val, digits = 3)),
                y = paste(pathway, "Score")
            ) +
            theme_classic()

        plot_list[[pathway]] <- p
    }
}

p_combined <- plot_grid(plotlist = plot_list, ncol = 2)
ggsave(file.path(results_dir, "Figures/04_KIF13A_Pathway_Correlations.pdf"), p_combined, width = 12, height = 10)

# ============================================================================
# 5. 生成机制简�?
# ============================================================================

sink(file.path(results_dir, "Reports/04_Mechanism_Summary.txt"))
cat("================================================\n")
cat("          KIF13A MECHANISM DISCOVERY            \n")
cat("================================================\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("1. Co-expression Network:\n")
cat("   - Identified Neighbors:", nrow(top_neighbors), "genes (R > 0.4)\n")
cat("   - Top 5 Neighbors:\n")
print(head(top_neighbors, 5))
cat("\n")

cat("2. Pathway Enrichment (Top 5 GO Terms):\n")
if (exists("ego") && !is.null(ego)) {
    print(head(as.data.frame(ego)[, c("ID", "Description", "p.adjust")], 5))
} else {
    cat("   No significant GO terms found.\n")
}
cat("\n")

cat("3. Pathway Correlations:\n")
for (pathway in names(plot_list)) {
    markers <- pathway_markers[[pathway]]
    valid_markers <- markers[markers %in% rownames(expr_log)]
    if (length(valid_markers) > 0) {
        pathway_score <- colMeans(expr_log[valid_markers, , drop = FALSE])
        cor_val <- cor(target_expr, pathway_score, method = "pearson")
        cat(paste("   -", pathway, ": R =", round(cor_val, 3), "\n"))
    }
}

sink()

message("Script 04 Completed. Mechanism exploration finished.")
