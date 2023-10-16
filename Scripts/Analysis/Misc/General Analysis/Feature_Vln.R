# Script information -------------------------------------------------------------------------------

# title: Seurat general plots - genes of interest
# author: James Favaloro
# date: 2021-09-08
# description: In this script we make some pretty pictures - Violin/Ridge/Feature/ plots


# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'Figures X' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scRepertoire)
  library(pheatmap)
})

# 2. Load data -------------------------------------------------------------------------------------

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("Data/black_list.rds")


# 3. Pre-processing --------------------------------------------------------------------------------

# Pull out the seurat object to analyse
x = integrated.list$all_ds.integrated

# Set the cluster resolution back to something more reasonable
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DefaultAssay(x) = "RNA"

# Check to ensure settings are correct
DimPlot(x)
VlnPlot(x, "CD69", split.by = "tissue")

# Backup x and remove the integrated list to free up memory
x_backup = x
rm(integrated.list)

# Split the seurat object into BM and PB

x_bm = subset(x, subset = tissue == "bm")
x_pb = subset(x, subset = tissue == "pb")

# Create directories if required
dir.create("Output/Figures/General/Vlnplots/")
dir.create("Output/Figures/General/Ridgeplots/")
dir.create("Output/Figures/General/Featureplots/")

# Set cluster colours to default seurat colours for clusters 0-6
# NB: It's good to set this as if a cluster isn't plotted colours get messed up
ccolours = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#619CFF", "#FF61C3")

# 4. User input ------------------------------------------------------------------------------------

# Store genes we're interested in plotting as a character vector we can iterate over
# Follow the convention: vector_of_genes = c("GENE1", "GENE2", "GENE3)

exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CD160", "CD244")
activation = c("AKT1", "AKT2", "AKT3", "ARAF", "BRAF", "CD247", "CD28", "CD3D", "CD3E", "CD3G", "CD80", "CD86", "CDC42", "CHUK", "CSK", "FOS",  "GRAP2", "HRAS", "IKBKB", "ITPR1", "JUN", "LAT", "LCK", "LCP2", "MAP2K1", "MAP2K2", "MAP3K1", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "NCK1", "NCK2", "NFKBIA", "NRAS", "PAK1", "PAK2", "PAK3", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PLCG1", "PPP3CA", "PPP3CB", "PPP3CC", "PRKCQ", "PTPRC", "RAC1", "RAF1", "SOS1", "SOS2",  "VAV1", "VAV2", "VAV3", "WAS", "ZAP70") #https://maayanlab.cloud/Harmonizome/gene_set/T+cell+activation/PANTHER+Pathways
effector = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "IFNG", "IL2")
regulators = c("EOMES", "TBX21", "TCF7", "BCL2", "TOX")
residency = c("CD69", "KLF2", "KLF3", "ITGA1", "S1PR1", "CD101", "ITGAE", "ICOS", "IRF4", "HOBIT")
homing = c("CCL3","CCL4","CCL5","CCL18","CCL23","CCL25","CCL28","CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCL12","CXCL13","CXCL16","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","XCL1","XCL2","CX3CR1")

#Set genes of interest (GOI) to whichever set of genes we want to plot
# i.e. if you want to plot effector genes: GOI = effector
GOI = homing

# Set Subset of interest (SOI) to whichever subset of cells we want to plot # Not yet implemented
SOI = x_dom

# 5.A Analysis - Violin plots ----------------------------------------------------------------------

# Violin plots
for (i in 1:length(GOI)) {
  temp = (VlnPlot(x, features = GOI[`i`], cols = ccolours, pt.size = 0))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Vlnplots/Vln_", GOI[i], ".png"))
}


# Violin plots, clusters split by tissue
for (i in 1:length(GOI)) {
  temp = (VlnPlot(x, features = GOI[`i`], split.by = "tissue", cols = c("red", "blue"), pt.size = 0))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Vlnplots/Tissue_", GOI[i], ".png"))
}

# Violin plots of BM clusters
for (i in 1:length(GOI)) {
  temp = (VlnPlot(x_bm, features = GOI[`i`], cols = ccolours, pt.size = 0))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Vlnplots/BM_", GOI[i], ".png"))
}

# Violin plots of PB clusters
for (i in 1:length(GOI)) {
  temp = (VlnPlot(x_pb, features = GOI[`i`], cols = ccolours, pt.size = 0))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Vlnplots/BM_", GOI[i], ".png"))
}

# 5.B Analysis - Ridge plots -----------------------------------------------------------------------

# Ridge plots
for (i in 1:length(GOI)) {
  temp = (RidgePlot(x, features = GOI[`i`], cols = ccolours))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Ridgeplots/Ridge_", GOI[i], ".png"))
}

# Ridge plots, split by tissue
for (i in 1:length(GOI)) {
  temp = (RidgePlot(x, features = GOI[`i`], group.by = "tissue", cols = c("red", "blue")))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Ridgeplots/Tissue_", GOI[i], ".png"))
}

# Ridge plots of BM clusters
for (i in 1:length(GOI)) {
  temp = (RidgePlot(x_bm, features = GOI[`i`], cols = ccolours))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Ridgeplots/BM_", GOI[i], ".png"))
}

# Ridge plots of PB clusters
for (i in 1:length(GOI)) {
  temp = (RidgePlot(x_pb, features = GOI[`i`], cols = ccolours))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Ridgeplots/PB_", GOI[i], ".png"))
}

# 5.C Analysis - Feature plots ---------------------------------------------------------------------

# Feature plots
for (i in 1:length(GOI)) {
  temp = (FeaturePlot(x, features = GOI[`i`], cols = c("lightgrey", "purple")))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Featureplots/Feature_", GOI[i], ".png"))
}

# Feature plots, split by tissue
for (i in 1:length(GOI)) {
  temp = (FeaturePlot(x, features = GOI[`i`], cols = c("lightgrey", "purple"), split.by = "tissue"))
  ggsave(width = 12.3, height = 7, units = "cm", plot = temp, filename = paste0("Output/Figures/General/Featureplots/Tissue_", GOI[i], ".png"))
}


# Feature plots of BM clusters
for (i in 1:length(GOI)) {
  temp = (FeaturePlot(x_bm, features = GOI[`i`], cols = c("lightgrey", "purple")))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Featureplots/BM_", GOI[i], ".png"))
}

# Feature plots of PB clusters
for (i in 1:length(GOI)) {
  temp = (FeaturePlot(x_pb, features = GOI[`i`], cols = c("lightgrey", "purple")))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Featureplots/PB_", GOI[i], ".png"))
}

# Overlay plots
FeaturePlot(x, features = c("GZMK", "PRF1"), split.by = "tissue", blend = TRUE)
FeaturePlot(x, features = c("GZMB", "PRF1"), split.by = "tissue", blend = TRUE)
FeaturePlot(x, features = c("GZMK", "CD69"), split.by = "tissue", blend = TRUE)
FeaturePlot(x, features = c("GZMB", "CD69"), split.by = "tissue", blend = TRUE)
#ggsave(dpi = "print", width = 12.3, height = 7, units = "cm", filename = ("Output/Figures/General/Featureplots/GZMK_PRF1.png"))
# Save as 1920x1080
