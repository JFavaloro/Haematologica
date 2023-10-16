# Script information -----------------------------------------------------------

# title: Analysis
# author: James Favaloro
# date: 2021-03-10
# description: In this script we make some pretty pictures - 2. Clusters heatmap

# THIS WILL BE SUPERCEEDED BY THE DE SCRIPTS

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
annotations = read.csv("Data/annotations.csv")

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

# Run FindAllMarkers to determine genes that define our clusters
markers = FindAllMarkers(x)

# Annotate the clusters using annotation sourced from WHERE
# NB - this will pull out genes we have excluded from clustering
annotated_markers <- inner_join(x = markers, 
                                y = annotations[, c("gene_name", "description")],
                                by = c("gene" = "gene_name")) %>% unique()

# Pull out the top 10 genes (based on level of expression) that define each cluster so we can visualise them

top_10_markers = annotated_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique = unique(top_10_markers$gene)

# Remove genes that don't exist in scaled data
top_10_unique_blacklisted = top_10_unique[top_10_unique %!in% black_list]

# Output figure 2.1 - Genes that define clusters in a heatmap
# NB: Package seemingly doesn't allow group.colours - have editied and saved
DoMultiBarHeatmap(x, features = top_10_unique_blacklisted, assay = "integrated", group.by = "seurat_clusters", additional.group.by = "tissue", group.bar = TRUE, group.colors = c("red", "blue"))

DoHeatmap(x, assay = "SCT", features = top_10_unique_blacklisted)

# Output figure 2.2 - Genes that define clusters 
DotPlot(x, assay = "SCT", features = top_10_unique_blacklisted, split.by = "tissue", cols = c("red", "blue")) + theme(axis.text.x = element_text(angle = 90))

