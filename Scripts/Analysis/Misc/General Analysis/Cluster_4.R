# Script information -------------------------------------------------------------------------------

# title: Cluster 4 deep dive
# author: James Favaloro
# date: 2021-09-10
# description: In this script we will interogate cluster 4 further.

# scRepertoire has changed colour palettes - need to resolve this

# https://www.frontiersin.org/articles/10.3389/fimmu.2018.01333/full

# Recent studies using MR1 tetramers revealed that not all of the TCR usage of human MAIT cells is restricted to TRAV1-2–TRAJ33
# Approximately, 30% of MR1-restricted cells are TRAV1-2 joined TRAJ20 or TRAJ12
# TRAV1-2–TRAJ33 are mostly paired with TRBV6-6 and TRBV20
# TRAV7-2–TRAJ33 is the dominant MAIT in PB
# TRAV7-2–TRAJ12 transcripts are higher in kidney and intestine biopsies from some individuals 

# MAIT cells also express PLZF
# Most MHC-restricted CD8+T cells express the CD8αβ heterodimer, but CD8+ MAIT cells express CD8αα homodimers, and some of them coexpress the CD8αβ heterodimer
# Human peripheral blood MAIT cells are CCR5+CCR6+CCR7−CCR9+/− CXCR3−CXCR4+/− CXCR6+

# Subset our data on the basis of gene expression and query cluster distribution
# NB: requires data.table package

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'Figures X' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scRepertoire)
  library(pheatmap)
  library(escape)
  library(dittoSeq)
  library(readxl)
  library(data.table)
})

# 2. Load data -------------------------------------------------------------------------------------

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("black_list.rds")

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))

# 3. Pre-processing --------------------------------------------------------------------------------

# pull out the seurat object to analyse
x = integrated.list$all_ds.integrated

# Set the cluster resolution back to something more reasonable
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DefaultAssay(x) = "RNA"

# Fix clonotype order - not sure why we need to do this but current code orders single and small in reverse(?)

slot(x, "meta.data")$cloneType <- factor(slot(x, "meta.data")$cloneType, 
                                         levels = c("Hyperexpanded (628 < X <= 6275)", 
                                                    "Large (63 < X <= 628)", 
                                                    "Medium (6 < X <= 63)", 
                                                    "Small (1 < X <= 6)", 
                                                    "Single (0 < X <= 1)", NA))

# Check to ensure settings are correct
DimPlot(x)
VlnPlot(x, "CD69", split.by = "tissue")

# Backup x and remove the integrated list to free up memory
x_backup = x
rm(integrated.list)

# Split the seurat object into BM and PB
x_bm = subset(x, subset = tissue == "bm")
x_pb = subset(x, subset = tissue == "pb")

x_4 = subset(x, seurat_clusters == "4")

vizVgenes()

vizGenes(df = x, gene = "V", chain = "TRA", separate = "tissue")

vizGenes(df = x_4, gene = "V", chain = "TRA", separate = "tissue")

# 4. Subset data -----------------------------------------------------------------------------------

# Classical/ Alternate MAIT
MAIT_1 = subset(x, subset = CTgene %like% "TRAV1-2.TRAJ33") # 95 cells!
MAIT_2 = subset(x, subset = CTgene %like% "TRAV1-2.TRAJ12") # 26 cells!
MAIT_3 = subset(x, subset = CTgene %like% "TRAV1-2.TRAJ20") # 360 cells!

# Come back to this after some reading - is 7-2 only in mice?
MAIT_4 = subset(x, subset = CTgene %like% "TRAV7-2.TRAJ12") # No cells!
MAIT_5 = subset(x, subset = CTgene %like% "TRAV7-2.TRAJ33") # No cells!

# Same with this
MAIT_6 = subset(x, subset = CTgene %like% "TRAV36.TRAJ34") # No cells!
MAIT_7 = subset(x, subset = CTgene %like% "TRAV36.TRAJ37") # No cells!

# NKT?
MAIT_8 = subset(x, subset = CTgene %like% "TRBV25") # No cells!

# Set subset of interest for downstream analysis

SOI = MAIT_3

# Split this into bm/pb

SOI_bm = subset(SOI, subset = tissue == "bm")
SOI_pb = subset(SOI, subset = tissue == "pb")

GOI = c("IL2", "IL3", "IL4", "IL13", "IL17", "PRF1", "GZMB", "CCL4", "VEGFA", "IFNG", "TNF", "CSF2")


# 5. Analyse ---------------------------------------------------------------------------------------

for (i in 1:length(GOI)) {
  temp = (VlnPlot(SOI, features = GOI[`i`], cols = ccolours, pt.size = 0))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Vlnplots/Vln_", GOI[i], ".png"))
}

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))

p1 = occupiedscRepertoire(x, x.axis = "cluster")
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1 + scale_color_manual(values = c("Hyperexpanded (628 < X <= 625)" = "#FF4B20", "Large (63 < X <= 628)" = "#FFB433", "Medium (6 < X <= 63)" = "C6FDEC", "Small (1 < X <= 6)" = "#7AC5FF", "Single (0 < X <=1)" = "0348A6"))

# BM distribution
p1 = occupiedscRepertoire(SOI_bm, x.axis = "cluster")
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1 + scale_fill_manual(values = colorblind_vector(5))

# PB distribution
p1 = occupiedscRepertoire(SOI_pb, x.axis = "cluster")
p1 <- p1 + scale_x_continuous(labels = c(0:6), breaks = c(0:6))
p1

y = expression2List(SOI, group = "cluster")
