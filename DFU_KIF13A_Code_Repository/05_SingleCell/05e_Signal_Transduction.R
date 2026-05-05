# ==============================================================================
# Script 47: Virtual Signal Transduction Modeling (In Silico Mechanosensing)
# Purpose: Mathematical Modeling of "The Deaf Cell" Hypothesis
# Logic:
#   Signal Strength = [Ligand]ext * [Receptor]surface * K_transport
#   Target: Demonstrate that KIF13A is the Rate-Limiting Step for GF Therapy
# Output: 3D Response Surface Plot (The "Dead Zone")
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
library(plotly)
library(viridis)

# Helper function for Seurat v5 compatibility
get_data <- function(obj) {
    if ("data" %in% Layers(obj)) {
        return(LayerData(obj, layer = "data"))
    } else {
        return(GetAssayData(obj, slot = "data"))
    }
}

# 2. Data Loading (Fast Track) --------------------------------------------
project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
analysis_dir <- file.path(project_dir, "所有真实数据分�?)
sc_dir <- file.path(analysis_dir, "Validation_Data/SingleCell")

message("Loading Single-Cell Data for Signal Modeling...")
files <- list.files(sc_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) files <- list.files(sc_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)

# Load Subset
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

# Preprocess
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & percent.mt < 15)
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, nfeatures = 2000)

message("Data loaded: ", ncol(sc_obj), " cells.")

# 3. Model Definition: The "Logistics Equation" ---------------------------
# We define Cellular Response (R) as a function of:
# L: External Ligand Concentration (e.g., PDGF-BB, VEGF) - Simulated Variable
# T: Transport Capacity (KIF13A level) - Measured Variable
# Rec: Receptor Expression (PDGFRB, VEGFR2) - Measured Variable

# Biological Logic:
# The amount of functional receptor on *Surface* = Total_Receptor * f(Transport)
# f(Transport) ~ Sigmoid(KIF13A)
# Response ~ Hill_Equation(Ligand) * Surface_Receptor

# Hill Function for Transport Efficiency
transport_efficiency <- function(kif13a_expr, k_half = 0.5, n = 2) {
    # Normalized KIF13A (0 to 1 range approx)
    x <- pmax(kif13a_expr, 0)
    # Hill equation: Cooperativity in transport
    return(x^n / (k_half^n + x^n))
}

# 4. Simulation on Real Data ---------------------------------------------
message("--- Running Signal Transduction Simulation ---")

# Step A: Extract Real Expression Values
# We focus on the PDGFRB pathway (classic healing signal)
# Receptor: PDGFRB
# Transporter: KIF13A
# Target Response Gene (Downstream): FOS (Immediate Early Gene)

genes_of_interest <- c("KIF13A", "PDGFRB", "KDR", "FOS", "JUN")
valid_genes <- genes_of_interest[genes_of_interest %in% rownames(sc_obj)]

# Robust Data Extraction
# Robust Data Extraction
expr_data <- tryCatch(
    {
        # Try V5 style first
        if (packageVersion("Seurat") >= "5.0.0") {
            if (!"data" %in% Layers(sc_obj)) sc_obj <- JoinLayers(sc_obj)
            as.matrix(GetAssayData(sc_obj, layer = "data")[valid_genes, , drop = FALSE])
        } else {
            as.matrix(GetAssayData(sc_obj, slot = "data")[valid_genes, , drop = FALSE])
        }
    },
    error = function(e) {
        # Fallback for older versions or if JoinLayers fails
        message("Fallback data access...")
        as.matrix(GetAssayData(sc_obj, layer = "data")[valid_genes, , drop = FALSE])
    }
)

# Step B: Create Simulation Grid
# We simulate a "Clinical Trial" of increasing Growth Factor doses
# Dose range: 0 to 100 (relative units)
doses <- seq(0, 5, length.out = 50)

# We simulate the response for EACH CELL at EACH DOSE
# To make it computable, we sample 500 representative cells (mix of Healed/Non-healed)
set.seed(42)
cells_to_sim <- sample(colnames(sc_obj), 500)
sim_data <- list()

for (cell in cells_to_sim) {
    kif <- expr_data["KIF13A", cell]
    # Receptor level (Strictly Real Data)
    if ("PDGFRB" %in% rownames(expr_data)) {
        rec <- expr_data["PDGFRB", cell]
    } else if ("KDR" %in% rownames(expr_data)) {
        rec <- expr_data["KDR", cell]
    } else {
        # If no key receptor is detected in this dataset, we cannot simulate signal for this cell
        rec <- NA
    }

    if (!is.na(rec)) {
        # Calculate Surface Receptor (The "Effective" Receptor)
        # If KIF13A is 0, Surface Receptor -> 0 (logistics failure)
        # We normalize KIF expression to 0-3 range for sigmoid based on observed max
        max_kif <- max(expr_data["KIF13A", ], na.rm = TRUE)
        if (max_kif == 0) max_kif <- 1 # Prevent division by zero
        kif_norm <- (kif / max_kif) * 3

        surf_rec <- rec * transport_efficiency(kif_norm)

        # Calculate Response to Dose D
        # Response = Vmax * Dose / (Km + Dose) * Surface_Receptor_Capacity
        # We assume 'Rec' defines the Vmax capacity
        responses <- (doses / (1 + doses)) * surf_rec

        sim_data[[cell]] <- data.frame(
            Cell = cell,
            KIF13A = kif,
            Receptor_Total = rec,
            Surface_Receptor = surf_rec,
            Dose = doses,
            Response = responses
        )
    }
}

sim_df <- do.call(rbind, sim_data)

# 5. Visualization: The "Therapeutic Dead Zone" --------------------------
# We aggregate to visualize the average response surface
# X: Dose (Treatment Intensity)
# Y: KIF13A Level (Transport Capacity)
# Z: Response (Healing Signal)

# Bin KIF13A into levels for surface plotting
sim_df$KIF_Bin <- cut(sim_df$KIF13A, breaks = 20)
surface_data <- sim_df %>%
    group_by(Dose, KIF_Bin) %>%
    summarise(
        Mean_Response = mean(Response),
        Mean_KIF = mean(KIF13A),
        .groups = "drop"
    )

# Prepare for ggplot heatmap / contour
pdf(file.path(analysis_dir, "Figures/47_Signal_Transduction_Model.pdf"), width = 10, height = 8)

# Plot A: The Response Surface (Heatmap Style)
p1 <- ggplot(surface_data, aes(x = Dose, y = Mean_KIF, fill = Mean_Response)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", name = "Healing Signal") +
    theme_minimal() +
    labs(
        title = "The 'Therapeutic Dead Zone'",
        subtitle = "Modeled Signal Transduction: Why Growth Factors Fail in KIF13A-Low Cells",
        x = "Exogenous Growth Factor Dose (Simulated Treatment)",
        y = "Intracellular KIF13A Level (Logistics Capacity)"
    ) +
    annotate("text", x = 3.5, y = 0.8, label = "Transport Threshold", color = "cyan", vjust = -0.5) +
    annotate("text", x = 4, y = 0.5, label = "Dead Zone:\nHigh Dose, No Response", color = "black", fontface = "bold") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "cyan")

# Plot B: Dose-Response Curves (Stratified)
# Split cells into High vs Low KIF13A
sim_df$Group <- ifelse(sim_df$KIF13A > median(sim_df$KIF13A), "High KIF13A (Transport-Competent)", "Low KIF13A (Transport-Deficient)")
curve_data <- sim_df %>%
    group_by(Dose, Group) %>%
    summarise(Mean_Resp = mean(Response), SE = sd(Response) / sqrt(n()))

p2 <- ggplot(curve_data, aes(x = Dose, y = Mean_Resp, color = Group)) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = Mean_Resp - SE, ymax = Mean_Resp + SE, fill = Group), alpha = 0.2, color = NA) +
    theme_bw() +
    scale_color_manual(values = c("dodgerblue", "firebrick")) +
    scale_fill_manual(values = c("dodgerblue", "firebrick")) +
    labs(
        title = "Predicted Treatment Efficacy",
        subtitle = "KIF13A Status Determines the Ceiling of Therapeutic Response",
        y = "Predicted Signal Output", x = "Growth Factor Dose (Arbitrary Units)"
    )

print(p1)
print(p2)
dev.off()

# Save Data
write.csv(surface_data, file.path(analysis_dir, "Tables/47_Signal_Model_Surface.csv"))

message("Signal Transduction Modeling Complete.")
message("Results saved to: Figures/47_Signal_Transduction_Model.pdf")
