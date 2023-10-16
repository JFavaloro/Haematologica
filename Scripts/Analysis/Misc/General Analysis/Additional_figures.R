# Script information -----------------------------------------------------------

# title: Addtional figures to finish this damn thesis off
# author: James Favaloro
# date: 2021-12-05
# description: In this script we make some pretty pictures

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'FINISH' script")

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

# Set Ident back to original
Idents(x) = "orig.ident"

# Plot Figure 1 QC

# Clean data - Number of genes detected
VlnPlot(x, "nCount_RNA", pt.size = 0, cols = c("red", "blue", "red", "blue")) + scale_y_continuous(limits = c(0,1E4))+ ggtitle("UMIs detected") + xlab("Sample") + ylab("UMIs detected") + NoLegend()
ggsave(filename = "Output/Figures/QC/Seurat_QC_nCount_RNA.pdf", device = "pdf", width = 10, height = 15, units = "cm")

# Clean data - Number of UMI detected
VlnPlot(x, "nFeature_RNA",  pt.size = 0, cols = c("red", "blue", "red", "blue")) + scale_y_continuous(limits = c(0,3000)) + ggtitle("Genes detected") + xlab("Sample") + ylab("Genes detected") + NoLegend()
ggsave(filename = "Output/Figures/QC/Seurat_QC_nFeature_RNA.pdf", device = "pdf", width = 10, height = 15, units = "cm")

# Clean data - Percentage Mitochondrial DNA
VlnPlot(x, "percent.mt",  pt.size = 0, cols = c("red", "blue", "red", "blue")) + scale_y_continuous(limits = c(0,20)) + ggtitle("% of mitochondrial genes") + xlab("Sample") + ylab("Mitochondrial gene %") + NoLegend()
ggsave(filename = "Output/Figures/QC/Seurat_QC_percent.mt.pdf", device = "pdf", width = 10, height = 15, units = "cm")

# Plot Figure 2 Find Variable Features

#NB: This is using pre-processed data - can come back to this if need be but just need a quick figure

# Use PB43 as an example
PB43 = all_seurat_filtered_SCT$pb43
PB43 = FindVariableFeatures(PB43, selection.method = "vst", nfeatures = 2000)
top20 = head(VariableFeatures(PB43), 20)
p1 = VariableFeaturePlot(PB43)
p2 = LabelPoints(plot = p1, points = top20, repel = TRUE)
p2
ggsave(filename = "Output/Figures/QC/Variable_features_pre_blacklist.png", device = "png")

# Blacklist genes and repeat figure - This doesn't seem to be working - use old figures for now...
PB43 = all_seurat_filtered_SCT$pb43
PB43 = FindVariableFeatures(PB43, selection.method = "vst", nfeatures = 2000)
top20 = head(top, 20)
p1 = VariableFeaturePlot(PB43)
p2 = LabelPoints(plot = p1, points = top20, repel = TRUE)
p2
ggsave(filename = "Output/Figures/QC/Variable_features_post_blacklist.png", device = "png")