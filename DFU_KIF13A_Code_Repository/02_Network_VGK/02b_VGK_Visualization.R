# ============================================================================
# DFU Real Data Analysis - Figure 1A Schematic (Concept Art)
# Script 02_Schematic: Generate High-Res "Virtual Knockout" Concept
# Target: Figure 1A (The "Hook" of the paper)
# ============================================================================

rm(list = ls())
gc()

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(igraph)
library(RColorBrewer)

project_dir <- "E:/人工智能学习/任海涛虚拟基因敲�?
results_dir <- file.path(project_dir, "所有真实数据分�?)
dir.create(file.path(results_dir, "Figures"), showWarnings = FALSE)

# ============================================================================
# Core Idea: Visualize "Hidden Hub" vs "Visible Leaves"
# ============================================================================
# Left: Standard View (Color = Fold Change). Hub is grey (low FC).
# Right: Topology View (Size = Impact). Hub is huge (High Impact).

pdf(file.path(results_dir, "Figures/01_VGK_Concept_Schematic.pdf"), width = 12, height = 6)

layout(matrix(1:2, 1, 2))
par(mar = c(2, 2, 4, 2))

# 1. Create a synthetic network
set.seed(123)
g <- make_star(30, mode = "undirected")
# Add some random connections to make it look biological
g <- add_edges(g, sample(2:30, 10, replace = T), sample(2:30, 10, replace = T))

V(g)$name <- c("KIF13A", paste0("G", 1:29))

# 2. Panel A-1: The "Blind Spot" of DEG (Traditional View)
# Scenario: Hub has low FC (Grey), Leaves have high FC (Red)
V(g)$color_deg <- c("grey80", sample(c("#E41A1C", "#377EB8"), 29, replace = T))
V(g)$size_deg <- 15 # All same size

plot(g,
    layout = layout_with_fr(g),
    vertex.color = V(g)$color_deg,
    vertex.size = V(g)$size_deg,
    vertex.label = NA,
    main = "Traditional DEG Analysis\n(Hub is Invisible)",
    sub = "Focus on Abundance (Color)"
)

legend("bottomleft",
    legend = c("High DE", "Low DE (Hub)"),
    col = c("#E41A1C", "grey80"), pch = 19, bty = "n", cex = 1.2
)

# 3. Panel A-2: The "Discovery" of VGK (Network View)
# Scenario: Hub is knocked out -> High Impact (Huge Size)
V(g)$color_vgk <- c("gold", rep("grey90", 29))
V(g)$size_vgk <- c(40, rep(10, 29)) # Hub is huge

plot(g,
    layout = layout_with_fr(g),
    vertex.color = V(g)$color_vgk,
    vertex.size = V(g)$size_vgk,
    vertex.label = ifelse(V(g)$name == "KIF13A", "KIF13A", NA),
    vertex.label.font = 2, vertex.label.cex = 1.2,
    main = "Virtual Gene Knockout (VGK)\n(Hub Revealed by Topology)",
    sub = "Focus on Impact (Size)"
)

legend("bottomleft",
    legend = c("Topological Driver", "Passenger"),
    col = c("gold", "grey90"), pch = 19, bty = "n", cex = 1.2
)

dev.off()

message("Figure 1A Schematic Generated.")
