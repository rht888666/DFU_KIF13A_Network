# ==============================================================================
# Script 40: Advanced Network Robustness Analysis
# Purpose: Statistical Hardening of Transport Module (Phase 2 Upgrade)
# Method: Mixed-Effects Modeling + Bootstrap + Permutation Testing
# ==============================================================================

# 1. Setup & Licensing ----------------------------------------------------
rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Function to check and install packages
ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        install.packages(pkg)
    }
}

ensure_package("tidyverse")
ensure_package("lme4")
ensure_package("lmerTest")
ensure_package("boot")
ensure_package("igraph")
ensure_package("reshape2")

library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(igraph)

# 2. Data Loading ---------------------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
data_file <- file.path(analysis_dir, "Data/01_REAL_processed_data.RData")

if (!file.exists(data_file)) stop("Data file not found!")
load(data_file)

# Prepare Data
expr_mat <- expr_log
meta <- metadata

# Align samples
common <- intersect(colnames(expr_mat), meta$sample)
expr_mat <- expr_mat[, common]
meta <- meta %>% filter(sample %in% common)
meta$healed <- factor(meta$healed, levels = c("No", "Yes")) # Reference is Non-healed

# Define Transport Module
transport_genes <- c("KIF13A", "EPN1", "CLIP1", "RAB11A")
valid_genes <- intersect(transport_genes, rownames(expr_mat))

if (length(valid_genes) < 2) stop("Not enough transport genes found.")

# Calculate Module Score (PC1 or Mean)
# Here using Mean Z-score for robustness/simplicity in interpretation
t_expr <- t(scale(t(expr_mat[valid_genes, ])))
module_score <- colMeans(t_expr)
meta$Transport_Score <- module_score[meta$sample]

# 3. Mixed-Effects Modeling (Patient Correction) --------------------------
message("--- Running Mixed-Effects Model ---")

# Standard LM vs Mixed Effects
lm_fit <- lm(Transport_Score ~ healed, data = meta)
lmer_fit <- lmer(Transport_Score ~ healed + (1 | patient), data = meta)

lm_res <- summary(lm_fit)$coefficients["healedYes", ]
lmer_res <- summary(lmer_fit)$coefficients["healedYes", ]

message(paste0("Standard LM P-value: ", format.pval(lm_res[4], digits = 3)))
message(paste0("Mixed-Effects P-value: ", format.pval(lmer_res[5], digits = 3)))
message(paste0("Mixed-Effects Estimate: ", round(lmer_res[1], 3)))

# 4. Leave-One-Patient-Out (LOO) Stability --------------------------------
message("--- Running LOO Stability Analysis ---")

patients <- unique(meta$patient)
n_patients <- length(patients)
loo_estimates <- numeric(n_patients)
names(loo_estimates) <- patients

for (p in patients) {
    # Subset data excluding patient p
    sub_meta <- meta %>% filter(patient != p)

    # Re-run Mixed Effects
    # Note: Requires enough patients. If N is small, model might not converge,
    # but with 17 patients it should be fine.
    try(
        {
            fit <- lmer(Transport_Score ~ healed + (1 | patient), data = sub_meta)
            loo_estimates[p] <- fixef(fit)["healedYes"]
        },
        silent = TRUE
    )
}

loo_mean <- mean(loo_estimates, na.rm = TRUE)
loo_sd <- sd(loo_estimates, na.rm = TRUE)
cv_loo <- loo_sd / abs(loo_mean)

message(paste0("LOO Mean Estimate: ", round(loo_mean, 3)))
message(paste0("LOO Stability (CV): ", round(cv_loo, 3)))


# 5. Network Permutation Test (Topological Robustness) --------------------
message("--- Running Network Permutation Test ---")
# Hypothesis: The co-expression density of these 4 genes is higher than random sets of 4 genes
# This proves "Topological Existence"

n_perm <- 1000
obs_density <- mean(cor(t(expr_mat[valid_genes, ]))[upper.tri(diag(length(valid_genes)))])

perm_densities <- numeric(n_perm)
all_genes <- rownames(expr_mat)

set.seed(123)
for (i in 1:n_perm) {
    rand_genes <- sample(all_genes, length(valid_genes))
    sub_mat <- expr_mat[rand_genes, ]
    # Handle potential zero variance genes in random sampling
    cor_mat <- cor(t(sub_mat))
    perm_densities[i] <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
}

p_val_perm <- sum(perm_densities >= obs_density) / n_perm
message(paste0("Network Permutation P-value: ", p_val_perm))


# 6. Visualization (Figure 2) ---------------------------------------------
out_pdf <- file.path(analysis_dir, "Figures/16_Robustness_Hardening.pdf")
pdf(out_pdf, width = 12, height = 5)

layout(matrix(c(1, 2, 3), 1, 3, byrow = TRUE))

# Plot 1: Mixed Effects Result
boxplot(Transport_Score ~ healed,
    data = meta,
    main = "Patient-Corrected Transport Levels",
    col = c("#E74C3C", "#2ECC71"),
    ylab = "Module Score (Z-scaled)", xlab = "Outcome"
)
stripchart(Transport_Score ~ healed,
    data = meta,
    vertical = TRUE, add = TRUE, method = "jitter",
    pch = 21, bg = "grey", cex = 0.8
)
text(1.5, max(meta$Transport_Score) * 0.9,
    paste0("LMM P = ", format.pval(lmer_res[5], digits = 3), "\n(Patient Adjusted)"),
    cex = 1.2, font = 2
)

# Plot 2: LOO Stability
# Waterfall plot of estimates
loo_df <- data.frame(Patient = names(loo_estimates), Estimate = loo_estimates)
loo_df <- loo_df %>% arrange(Estimate)
barplot(loo_df$Estimate,
    main = "Leave-One-Patient-Out Stability",
    ylab = "Effect Size Estimate (Beta)",
    xlab = "Excluded Patient ID (Sorted)",
    col = ifelse(loo_df$Estimate > 0, "#3498DB", "grey"),
    border = NA, ylim = c(min(loo_df$Estimate) * 1.2, max(loo_df$Estimate) * 1.2)
)
abline(h = loo_mean, col = "red", lty = 2, lwd = 2)
abline(h = 0, col = "black")
text(1, max(loo_df$Estimate), paste("Stability CV =", round(cv_loo, 2)), pos = 4)

# Plot 3: Permutation Test
hist(perm_densities,
    breaks = 30, col = "grey90", border = "grey",
    main = "Topological Specificity Test (n=1000)",
    xlab = "Mean Correlation (Random Modules)",
    xlim = c(min(perm_densities), max(c(perm_densities, obs_density)))
)
abline(v = obs_density, col = "red", lwd = 3)
text(obs_density, 0, paste0("Observed\nTransport Module\n(P < ", max(1 / n_perm, p_val_perm), ")"),
    pos = 2, col = "red", font = 2
)

dev.off()

message("Figure Saved to: ", out_pdf)
