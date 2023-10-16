# Script information -----------------------------------------------------------

# title: Analysis
# author: James Favaloro
# date: 2021-03-10
# description: In this script we make some pretty pictures - 3. scRepertoire

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

# 5.1 Have a look at clonotypic distribution

# Output figure 3.1 - expansion of clonotypes
# Have used NoLegend for time being until I can figure out how to customise the labels
DimPlot(x, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend() + ggtitle("Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot.png"), device = "png")

# Output figure 3.2 - expansion of clonotypes by tissue 
#NB Need no hyperexpanded clones in the BM - need manual colours
DimPlot(x_bm, group.by = "cloneType") +
  scale_color_manual(values = c("#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"), na.value="grey") + NoLegend() + ggtitle("BM Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_bm.png"), device = "png")

DimPlot(x_pb, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend() + ggtitle("PB Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_pb.png"), device = "png")

# Output figure 3.3 - expansion of clonotypes by sample
temp = subset(x, subset = orig.ident == "bm43")
DimPlot(temp, group.by = "cloneType") +
  scale_color_manual(values = c("#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"), na.value="grey") + NoLegend() + ggtitle("BM 43 Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_bm43.png"), device = "png")

temp = subset(x, subset = orig.ident == "pb43")
DimPlot(temp, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend() + ggtitle("PB 43 Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_pb43.png"), device = "png")

temp = subset(x, subset = orig.ident == "bm63")
DimPlot(temp, group.by = "cloneType") +
  scale_color_manual(values = c("#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"), na.value="grey") + NoLegend() + ggtitle("BM 63 Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_bm63.png"), device = "png")

temp = subset(x, subset = orig.ident == "pb63")
DimPlot(temp, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend() + ggtitle("PB 63 Clonotype Distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Dim_plot_pb63.png"), device = "png")

# Visualise distrubtion of clonotypes
p1 = occupiedscRepertoire(x, x.axis = "cluster", label = FALSE)
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1 + scale_fill_manual(values = colorblind_vector(5)) + NoLegend() + ggtitle("Clonotype Distribution by Cluster") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Barplot_distribution.png"), device = "png")

# BM distribution # Ask Sam RE: colours
p1 = occupiedscRepertoire(x_bm, x.axis = "cluster", label = FALSE)
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1 + scale_fill_manual(values = colorblind_vector(5)) + NoLegend() + ggtitle("BM Clonotype Distribution by Cluster") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Barplot_bm_distribution.png"), device = "png")

# PB distribution
p1 = occupiedscRepertoire(x_pb, x.axis = "cluster", label = FALSE)
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1 + scale_fill_manual(values = colorblind_vector(5)) + NoLegend() + ggtitle("PB Clonotype Distribution by Cluster") + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Barplot_pb_distribution.png"), device = "png")

# Per patient distribution # Ask Sam - why is NA not showing
occupiedscRepertoire(x, x.axis = "orig.ident", na.include = TRUE) + scale_fill_manual(values = colorblind_vector(5), na.value = "grey") + NoLegend()
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Barplot_sample_distribution.png"), device = "png")


# Diversity metrics for plotting
combined2 <- expression2List(x, split.by = "cluster")
clonalDiversity(combined2, exportTable = T)

# Are BM clusters more diverse than PB clusters?
combined2 <- expression2List(x_bm, split.by = "cluster")
clonalDiversity(combined2, exportTable = T)

# Are PB clusters more diverse than BM clusters?
combined2 <- expression2List(x_pb, split.by = "cluster")
clonalDiversity(combined2, exportTable = T)

# Overlap analysis by cluster
combined2 <- expression2List(x, split.by = "cluster")
clonalOverlap(combined2, method = "morisita")
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Morisita.png"), device = "png")

# BM overlap
combined2 <- expression2List(x_bm, split.by = "cluster")
clonalOverlap(combined2, method = "morisita")
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Morisita_BM.png"), device = "png")

# PB overlap
combined2 <- expression2List(x_pb, split.by = "cluster")
clonalOverlap(combined2, method = "morisita")
ggsave(filename = paste0("Output/Figures/General/scRepertoire/Morisita_PB.png"), device = "png")


# Merged overlap

#Rename Idents
x_bm = RenameIdents(x_bm, '0' = "0_bm", '1' = "1_bm", "2" = "2_bm", '3' = "3_bm", '4' = "4_bm", '5' = "5_bm", '6' = "6_bm")

x_pb = RenameIdents(x_pb, '0' = "0_pb", '1' = "1_pb", "2" = "2_pb", '3' = "3_pb", '4' = "4_pb", '5' = "5_pb", '6' = "6_pb")

# Merge Seurat object
x_merged = merge(x_bm, x_pb)

# Copy reductions and graphs?
x_merged@graphs = x@graphs
x_merged@reductions = x@reductions
x_merged@commands = x@commands
x_merged@tools = x@tools

# Run clonal overlap
Clonalnetwork = clonalNetwork(x_merged, reduction = "umap", identity = "cluster", cloneCall = "aa", exportTable = TRUE)
Clonalnetwork_figure = clonalNetwork(x_merged, reduction = "umap", identity = "cluster", cloneCall = "aa")
