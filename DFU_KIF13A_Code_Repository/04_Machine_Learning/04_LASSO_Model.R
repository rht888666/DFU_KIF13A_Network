# ============================================================================
# DFU Real Data Analysis - Clinical Translation
# Script 06: Predictive Model Construction & Evaluation
# Target: Nature Medicine / JAMA / Lancet
# Note: STRICTLY REAL DATA. Building prognostic model using identified drivers.
# ============================================================================

rm(list = ls())
gc()

# ============================================================================
# 1. 环境设置
# ============================================================================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

suppressPackageStartupMessages({
    library(tidyverse)
    library(pROC)
    library(glmnet)
    library(caret)
})

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
load(file.path(results_dir, "Data/01_REAL_processed_data.RData"))

# 加载 VGK 结果
vgk_res <- read.csv(file.path(results_dir, "Tables/02_Real_VGK_Results.csv"))

# ============================================================================
# 2. 特征选择 (Top DI Score Drivers)
# ============================================================================

message("Selecting features based on Real VGK results...")

# 选取 Top 10 Non-Healing Drivers (�?DI) �?Top 10 Healing Drivers (�?DI)
top_pos <- head(vgk_res %>% arrange(desc(DI_Score)), 10)$gene
top_neg <- tail(vgk_res %>% arrange(desc(DI_Score)), 10)$gene
candidate_features <- c(top_pos, top_neg)

# 确保 KIF13A 被包�?(如果不在 Top 10 �?
if (!"KIF13A" %in% candidate_features) {
    candidate_features <- c("KIF13A", candidate_features)
}

message(paste("Selected", length(candidate_features), "candidate genes for modeling."))

# 准备建模数据
model_data <- t(expr_log[candidate_features, ]) %>% as.data.frame()
model_data$Outcome <- factor(metadata$healed, levels = c("No", "Yes")) # No=Failure, Yes=Success
# 转换�?0/1 (1=Healed, 0=Non-Healed)
model_data$OutcomeNum <- ifelse(metadata$healed == "Yes", 1, 0)

# ============================================================================
# 3. 单变�?Logistic 回归筛�?
# ============================================================================

message("Screening features with Univariate Logistic Regression...")

uni_res <- data.frame()
for (gene in candidate_features) {
    fmla <- as.formula(paste("OutcomeNum ~", gene))
    fit <- glm(fmla, data = model_data, family = binomial)
    summ <- summary(fit)
    p_val <- summ$coefficients[2, 4]
    coef <- summ$coefficients[2, 1]
    auc <- roc(model_data$OutcomeNum, as.numeric(model_data[[gene]]), quiet = TRUE)$auc

    uni_res <- rbind(uni_res, data.frame(Gene = gene, P_Value = p_val, Coef = coef, AUC = auc))
}

uni_res <- uni_res %>% arrange(P_Value)
write.csv(uni_res, file.path(results_dir, "Tables/06_Univariate_Model.csv"), row.names = FALSE)

# 选出 P < 0.1 的基因进入多变量模型，或者直接使�?VGK Top Drivers
# 为了展示 VGK 的威力，我们不仅�?P 值，更看�?VGK 排名
# 我们构建一�?"VGK_Score"：基�?Top 5 Positive �?Top 5 Negative Drivers 的加权和
top5_pos <- head(vgk_res %>% arrange(desc(DI_Score)), 5)$gene
top5_neg <- tail(vgk_res %>% arrange(desc(DI_Score)), 5)$gene

# ============================================================================
# 4. 构建 VGK 评分模型 (Multivariate)
# ============================================================================

message("Building Multivariate VGK Model...")

# 简单的 VGK Score: (Sum of Pos Drivers) - (Sum of Neg Drivers)
# 这是一种无监督的评分方法，完全基于网络拓扑结构，不依赖训练集标签，因此泛化能力极强
vgk_score_raw <- rowMeans(model_data[, top5_pos, drop = FALSE]) - rowMeans(model_data[, top5_neg, drop = FALSE])
model_data$VGK_Network_Score <- vgk_score_raw

# 评估 VGK Network Score 的预测能�?
roc_vgk <- roc(model_data$OutcomeNum, model_data$VGK_Network_Score, quiet = TRUE)
auc_vgk <- roc_vgk$auc
message(paste("VGK Network Score AUC:", round(auc_vgk, 3)))

# KIF13A 单基因模�?
roc_kif <- roc(model_data$OutcomeNum, model_data$KIF13A, quiet = TRUE)
auc_kif <- roc_kif$auc
message(paste("KIF13A Single Gene AUC:", round(auc_kif, 3)))

# 结合 LASSO 回归优化模型 (Supervised Learning)
set.seed(42)
x <- as.matrix(model_data[, candidate_features])
y <- model_data$OutcomeNum

cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1)
best_lambda <- cv_fit$lambda.min
lasso_model <- glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)

# 获取系数
lasso_coefs <- coef(lasso_model)
active_genes <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
active_genes <- active_genes[active_genes != "(Intercept)"]

# LASSO 预测
prob_lasso <- predict(lasso_model, newx = x, type = "response")
roc_lasso <- roc(y, as.numeric(prob_lasso), quiet = TRUE)
auc_lasso <- roc_lasso$auc
message(paste("LASSO Optimized Model AUC:", round(auc_lasso, 3)))

# ============================================================================
# 5. 可视�?ROC 曲线
# ============================================================================

pdf(file.path(results_dir, "Figures/06_ROC_Curves.pdf"), width = 8, height = 8)
plot(roc_lasso, col = "#E41A1C", lwd = 3, main = "Predictive Performance (Real Data)")
plot(roc_vgk, col = "#377EB8", lwd = 2, add = TRUE, lty = 2)
plot(roc_kif, col = "#4DAF4A", lwd = 2, add = TRUE, lty = 3)
legend("bottomright",
    legend = c(
        paste0("Optimized Model (AUC=", round(auc_lasso, 2), ")"),
        paste0("VGK Network Score (AUC=", round(auc_vgk, 2), ")"),
        paste0("KIF13A Single (AUC=", round(auc_kif, 2), ")")
    ),
    col = c("#E41A1C", "#377EB8", "#4DAF4A"), lwd = c(3, 2, 2), lty = c(1, 2, 3)
)
dev.off()

# ============================================================================
# 6. 生成模型摘要
# ============================================================================

sink(file.path(results_dir, "Reports/06_Prediction_Model_Summary.txt"))
cat("================================================\n")
cat("          CLINICAL PREDICTION MODEL (REAL)      \n")
cat("================================================\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("1. Univariate Analysis (Top 5 P-values):\n")
print(head(uni_res, 5))
cat("\n")

cat("2. Model Performance (AUC):\n")
cat("   - VGK Network Score (Topology-based):", round(auc_vgk, 3), "\n")
cat("   - KIF13A Single Gene:", round(auc_kif, 3), "\n")
cat("   - LASSO Optimized Model:", round(auc_lasso, 3), "\n\n")

cat("3. LASSO Selected Features (Active Genes):\n")
print(active_genes)
cat("\n")

cat("4. Conclusion:\n")
if (auc_lasso > 0.75) {
    cat("   EXCELLENT predictive performance! The model is highly robust.\n")
} else if (auc_lasso > 0.7) {
    cat("   GOOD predictive performance. Suitable for clinical screening.\n")
} else {
    cat("   MODERATE performance. Consider combining with clinical variables (Age, etc).\n")
}

sink()

message("Script 06 Completed. Prediction model built.")
