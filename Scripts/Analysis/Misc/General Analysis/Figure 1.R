# Script information -------------------------------------------------------------------------------

# title: Analysis
# author: James Favaloro
# date: 2021-03-10
# description: In this script we make some pretty pictures - 1. Clusters

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'Figure 1' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(scRepertoire)
  library(pheatmap)
  library(escape)
  library(dittoSeq)
  library(DoMultiBarHeatmap)
})

# 2. Load data, define functions and colours -------------------------------------------------------

# Load data

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("Data/black_list.rds")
annotation = read.csv("Data/annotation.csv")

# Define functions

# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Define colours

# Define colour scheme for scRepertoire
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
# Define colour scheme for Seurat
seurat_colours = c("0" = "#F8766D", "1" = "#D39200", "2" = "#93AA00", "3" = "#00BA38", "4" = "#00C19F", "5" = "#619CFF", "6" = "#FF61C3")

# 3. Pre-processing  -------------------------------------------------------------------------------

# Pull out the seurat object to analyse
x = integrated.list$all_ds.integrated

# Set the cluster resolution back to something more reasonable
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DefaultAssay(x) = "RNA"

# Fix clonotype order - not sure why but current code orders single and small in reverse(?)
slot(x, "meta.data")$cloneType <- factor(slot(x, "meta.data")$cloneType, 
                                         levels = c("Hyperexpanded (628 < X <= 6275)", 
                                                    "Large (63 < X <= 628)", 
                                                    "Medium (6 < X <= 63)", 
                                                    "Small (1 < X <= 6)", 
                                                    "Single (0 < X <= 1)", NA))

# Backup x and remove the integrated list to free up memory
x_backup = x
rm(integrated.list)

# Split the seurat object into BM and PB
x_bm = subset(x, subset = tissue == "bm")
x_pb = subset(x, subset = tissue == "pb")

# Create directories if required
dir.create("Output/Figures/General/Seurat/")
dir.create("Output/Figures/General/scRepertoire/")
dir.create("Output/Figures/General/Vlnplots/")
dir.create("Output/Figures/General/Ridgeplots/")
dir.create("Output/Figures/General/Featureplots/")


# 4. User input ------------------------------------------------------------------------------------

# Store genes we're interested in plotting as a character vector we can iterate over
# Follow the convention: vector_of_genes = c("GENE1", "GENE2", "GENE3)

exhaustion = c("PDCD1", "HAVCR2", "LAG3", "CD244", "ENTPD1", "TIGIT", "CD160", "BTLA")
activation = c("GZMA")
effector = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "EOMES", "TBET")
regulators = c("EOMES", "TBET", "TCF7", )
residency = c("CD69", )

#Set genes of interest (GOI) to whichever set of genes we want to plot
# i.e. if you want to plot effector genes: GOI = effector
GOI = effector

# Set Subset of interest (SOI) to whichever subset of cells we want to plot # Not yet implemented
SOI = x


# 5. Analyse  --------------------------------------------------------------------------------------

# 4.1 Look at cluster disribution
# Output images for figure 1 - cluster distribution and quantification

# Output figure 1.1 - clusters
DimPlot(x) + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot.png"), device = "png")

# Output figure 1.2 - clusters split by/ grouped by tissue
DimPlot(x_bm) + ggtitle("BM Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_bm.png"), device = "png")
DimPlot(x_pb) + ggtitle("PB Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_pb.png"), device = "png")
DimPlot(x, group.by = "tissue", cols = c("red", "blue")) + ggtitle("Tissue Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_tissue.png"), device = "png")

# Output figure 1.3 - clusters split by/ grouped by sample
temp = subset(x, subset = orig.ident == "bm43")
DimPlot(temp) + ggtitle("BM 43 Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_bm43.png"), device = "png")
temp = subset(x, subset = orig.ident == "pb43")
DimPlot(temp) + ggtitle("PB 43 Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_pb43.png"), device = "png")
temp = subset(x, subset = orig.ident == "bm63")
DimPlot(temp) + ggtitle("BM 63 Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_bm63.png"), device = "png")
temp = subset(x, subset = orig.ident == "pb63")
DimPlot(temp) + ggtitle("PB 63 Clusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/Seurat/Dim_plot_pb63.png"), device = "png")

# Output figure 1.4 - cluster quantification split by tissue as heatmap/barplot # Ask SAM for help
y = table(x@active.ident, x@meta.data$tissue)
png("Output/Figures/General/Seurat/Heatmap_tissue.png", width = 1500, height = 1500, units = "px")
pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20, number_format = "%.0f")
dev.off()
png("Output/Figures/General/Seurat/Barplot_tissue.png")
barplot(y, col = seurat_colours, legend.text = FALSE, yaxp=c(0,12500, 4))
dev.off()

# Output figure 1.5 - cluster quantification split by sample as heatmap/barplot
y = table(x@active.ident, x@meta.data$orig.ident)
pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20, number_format = "%.0f")
barplot(y, col = seurat_colours, legend.text = FALSE, yaxp=c(0,6250, 4))
