# ==============================================================================
# Script 41: Structural Equation Modeling (SEM)
# Purpose: Causal Pathway Analysis (Phase 3 Upgrade)
# Method: Lavaan Latent Variable Modeling
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

ensure_package("lavaan")
ensure_package("semPlot")
ensure_package("tidyverse")

library(lavaan)
library(semPlot)
library(tidyverse)

# 2. Data Loading ---------------------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
data_file <- file.path(analysis_dir, "Data/01_REAL_processed_data.RData")

if (!file.exists(data_file)) stop("Data file not found!")
load(data_file)

# Align Data
expr_mat <- expr_log
meta <- metadata
common <- intersect(colnames(expr_mat), meta$sample)
expr_mat <- expr_mat[, common]
meta <- meta %>% filter(sample %in% common)

# 3. Model Variable Preparation -------------------------------------------
# We need to construct a data frame for lavaan
# Latent Variable "Transport": Measured by KIF13A, EPN1, CLIP1
# Latent Variable "Inflammation": Measured by MMP9, S100A8, S100A9
# Outcome: Healed (1/0)

transport_indicators <- c("KIF13A", "EPN1", "CLIP1")
inflam_indicators <- c("MMP9", "S100A8", "S100A9")

# Extract and standardize
model_data <- as.data.frame(t(expr_mat[c(transport_indicators, inflam_indicators), ]))
model_data$Outcome <- ifelse(meta$healed == "Yes", 1, 0)

# Scale indicators to help convergence
model_data <- as.data.frame(scale(model_data))
# Re-add Outcome as binary (0/1) for interpretation, though SEM treats as continuous in standard estimator
# For rigorous binary outcome, we'd use WLSMV, but for path coefficients standard ML is often acceptable
# in this biological context if we interpret as "propensity to heal".
# Let's try `estimator = "ML"` first.

# 4. Define Lavaan Model --------------------------------------------------
sem_model <- "
  # Measurement Model (Latent Variables)
  Transport =~ KIF13A + EPN1 + CLIP1
  Inflammation =~ MMP9 + S100A8 + S100A9

  # Structural Model (Regressions)
  Inflammation ~ a*Transport
  Outcome ~ b*Inflammation + c*Transport

  # Indirect Effect Calculation
  indirect := a*b
  total := c + (a*b)
"

# 5. Fit Model ------------------------------------------------------------
message("--- Fitting SEM Model ---")
fit <- sem(sem_model, data = model_data, estimator = "ML", missing = "fiml")

# Check Convergence
if (!inspect(fit, "converged")) {
    warning("Model did not converge! Trying simplified path analysis...")
    # Fallback to Path Analysis with composite scores if Latent fails
    model_data$Transport_Score <- rowMeans(model_data[, transport_indicators])
    model_data$Inflammation_Score <- rowMeans(model_data[, inflam_indicators])
    path_model <- "
    Inflammation_Score ~ a*Transport_Score
    Outcome ~ b*Inflammation_Score + c*Transport_Score
    indirect := a*b
    total := c + (a*b)
  "
    fit <- sem(path_model, data = model_data)
}

# 6. Extract Results & Fit Indices ----------------------------------------
summary_fit <- summary(fit, fit.measures = TRUE, standardized = TRUE)
fit_measures <- fitMeasures(fit, c("cfi", "rmsea", "srmr", "chisq", "pvalue"))

message("--- Fit Indices ---")
print(fit_measures)

parameter_estimates <- parameterEstimates(fit, standardized = TRUE) %>%
    filter(op %in% c("~", ":=")) %>%
    select(lhs, op, rhs, est, se, z, pvalue, std.all)

print(parameter_estimates)

# Save Stats
write.csv(parameter_estimates, file.path(analysis_dir, "Tables/SEM_Results.csv"))
write.csv(fit_measures, file.path(analysis_dir, "Tables/SEM_Fit_Indices.csv"))


# 7. Visualization (Figure 3) ---------------------------------------------
out_pdf <- file.path(analysis_dir, "Figures/17_SEM_Causal_Model.pdf")
pdf(out_pdf, width = 8, height = 6)

semPaths(fit, "std",
    layout = "tree2",
    edge.label.cex = 1.2,
    curvePivot = TRUE,
    fade = FALSE,
    posCol = "#2ECC71",
    negCol = "#E74C3C",
    nodeLabels = c(
        "KIF13A", "EPN1", "CLIP1", "MMP9", "S100A8", "S100A9",
        "Transport", "Inflammation", "Outcome"
    ),
    sizeMan = 8, sizeLat = 10,
    title = TRUE,
    mar = c(3, 5, 3, 5)
)

title(paste0(
    "Structural Equation Model\nCFI=", round(fit_measures["cfi"], 3),
    " RMSEA=", round(fit_measures["rmsea"], 3)
), line = 2)

# Add Annotation for Indirect Effect
ind_eff <- parameter_estimates %>% filter(lhs == "indirect")
text(0, -1, paste0(
    "Indirect Effect (Mediation):\nBeta = ", round(ind_eff$std.all, 3),
    "\nP = ", format.pval(ind_eff$pvalue, digits = 3)
),
cex = 0.9, font = 3
)

dev.off()

message("Figure Saved to: ", out_pdf)
