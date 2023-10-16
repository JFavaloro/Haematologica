# Script information -------------------------------------------------------------------------------

# Title: 10x dataset - DE Analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will import the previous analysis, annotate the clusters we have
# identified and create some additional figures. We will then export the Seurat objects for further
# downstream analysis. Analysis has revealed the 'correct' number of clusters is 6 for the all_10x
# and BM_10x object and 5 clusters for the PB_10x object. We will standardise the colours for
# clusters and ensure that plots output in a sensible manner.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. 10x dataset - DE Analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, openxlsx, EnhancedVolcano, scales, DoMultiBarHeatmap)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data
integrated_list <- readRDS("Data/R_out/10x/integrated_list_10_PC.rds")

# Load previous analysis
load("Data/R_out/10x/DE_Workup_10x.rds")

# Define cluster names for naming of lists
cluster_names_full <- 
  c("Combined ~ T[EM]", "Combined ~ T[TE]", "Combined ~ T[N]", 
    "Combined ~ Cyto ~ T[EM]", "Combined ~ P[RE] ~ EX", "Combined ~ T[CM]", 
    "BM ~ T[EM]", "BM ~ T[N]", "BM ~ T[TE]", 
    "BM ~ T[CM]", "BM ~ P[RE] ~ EX", "BM ~ Cyto ~ T[EM]", 
    "PB ~ T[TE]", "PB ~ T[EM]", "PB ~ T[N]", 
    "PB ~ T[CM]", "PB ~ Cyto ~ T[EM]")

excel_cluster_names_full <- 
  c("all_10x_TEM", "all_10x_TTE", "all_10x_TN", 
    "all_10x_Cyto_TEM", "all_10x_PRE_EX", "all_10x_TCM", 
    "BM_10x_TEM", "BM_10x_TN", "BM_10x_TTE", 
    "BM_10x_TCM", "BM_10x_PRE_EX", "BM_10x_Cyto_TEM", 
    "PB_10x_TTE", "PB_10x_TEM", "PB_10x_TN", 
    "PB_10x_TCM", "PB_10x_Cyto_TEM")

cluster_names_all_10x <- c("T[EM]", "T[TE]","T[N]", "Cyto ~ T[EM]", "P[RE] ~ EX", "T[CM]")
cluster_names_BM_10x <- c("T[EM]", "T[N]", "T[TE]", "T[CM]", "P[RE] ~ EX", "Cyto ~ T[EM]")
cluster_names_PB_10x <- c("T[TE]", "T[EM]", "T[CM]", "T[N]", "Cyto ~ T[EM]")

excel_cluster_names_all_10x <- c("TEM", "TTE","TN", "Cyto_TEM", "PRE_EX", "TCM")
exceL_cluster_names_BM_10x <- c("TEM", "TN", "TTE", "TCM", "PRE_EX", "Cyto_TEM")
excel_cluster_names_PB_10x <- c("TTE", "TEM", "TCM", "TN", "Cyto_TEM")

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")

scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))

# Create character vector for plot colours
all_10x_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[TE]" = "#B79F00", "T[N]" = "#00BA38", 
    "Cyto ~ T[EM]" = "#00BFC4", "P[RE] ~ EX" = "#619CFF", "T[CM]" = "#F564E3")

BM_10x_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[N]" = "#00BA38", "T[TE]" = "#B79F00",
    "T[CM]" = "#F564E3", "P[RE] ~ EX" = "#619CFF", "Cyto ~ T[EM]" = "#00BFC4")

PB_10x_plot_cols <- 
  c("T[TE]" = "#B79F00", "T[EM]" = "#F8766D", "T[CM]" = "#F564E3", 
    "T[N]" = "#00BA38", "Cyto ~ T[EM]" = "#00BFC4")

ccolours_all_10x <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
ccolours_BM_10x <- c("#00BA38", "#F8766D", "#B79F00", "#F564E3", "#619CFF", "#00BFC4")
ccolours_PB_10x <- c("#B79F00", "#F8766D", "#F564E3", "#00BA38", "#00BFC4")


# 3. Pre-processing --------------------------------------------------------------------------------

# Fix names for "cloneType" to allow easier manipulation of data
integrated_list$all_ds.integrated@meta.data <- 
  integrated_list$all_ds.integrated@meta.data %>% 
  separate(cloneType, into = "cloneType", sep = "\\(")  
integrated_list$BM_ds.integrated@meta.data <- 
  integrated_list$BM_ds.integrated@meta.data %>% 
  separate(cloneType, into = "cloneType", sep = "\\(")  
integrated_list$PB_ds.integrated@meta.data <- 
  integrated_list$PB_ds.integrated@meta.data %>% 
  separate(cloneType, into = "cloneType", sep = "\\(")  

# Return character to factor for cloneType
integrated_list$all_ds.integrated@meta.data$cloneType <- 
  as.factor(integrated_list$all_ds.integrated@meta.data$cloneType)
integrated_list$BM_ds.integrated@meta.data$cloneType <- 
  as.factor(integrated_list$BM_ds.integrated@meta.data$cloneType)
integrated_list$PB_ds.integrated@meta.data$cloneType <- 
  as.factor(integrated_list$PB_ds.integrated@meta.data$cloneType)

# Extract the integrated objects we're interested in and set the correct number of clusters
all_10x <- integrated_list$all_ds.integrated
BM_10x <- integrated_list$BM_ds.integrated
PB_10x <- integrated_list$PB_ds.integrated

# Set resolution to 0.18 for all_10x object
all_10x$seurat_clusters <- 
  all_10x$integrated_snn_res.0.18
Idents(all_10x) <- "seurat_clusters"

# Set resolution to 0.28 for BM_10x object
BM_10x$seurat_clusters <- 
  BM_10x$integrated_snn_res.0.28
Idents(BM_10x) <- "seurat_clusters"

# Set resolution to 0.16 for PB_10x object
PB_10x$seurat_clusters <- PB_10x$integrated_snn_res.0.16
Idents(PB_10x) <- "seurat_clusters"

# 3.1 Rename the Integrated all_10x object ---------------------------------------------------------

all_10x$Cluster <- all_10x$seurat_clusters
Idents(all_10x) <- "Cluster"
names(cluster_names_all_10x) <- levels(all_10x)
all_10x <- RenameIdents(all_10x, cluster_names_all_10x)
levels(all_10x$Cluster) <- cluster_names_all_10x

# 3.2 Rename the Integrated BM_10x object ----------------------------------------------------------

BM_10x$Cluster <- BM_10x$seurat_clusters
Idents(BM_10x) <- "Cluster"
names(cluster_names_BM_10x) <- levels(BM_10x)
BM_10x <- RenameIdents(BM_10x, cluster_names_BM_10x)
levels(BM_10x$Cluster) <- cluster_names_BM_10x

# 3.3 Rename the Integrated PB_10x object ----------------------------------------------------------

PB_10x$Cluster <- PB_10x$seurat_clusters
Idents(PB_10x) <- "Cluster"
names(cluster_names_PB_10x) <- levels(PB_10x)
PB_10x <- RenameIdents(PB_10x, cluster_names_PB_10x)
levels(PB_10x$Cluster) <- cluster_names_PB_10x

# 3.4 Continue Pre-processing ----------------------------------------------------------------------

# Split the all_10x object into BM and PB objects
all_10x_BM <- subset(all_10x, subset = Tissue == "BM")
all_10x_PB <- subset(all_10x, subset = Tissue == "PB")

# Move all objects into a list
all_list <- list("all_10x" = all_10x, "BM_10x" = BM_10x, "PB_10x" = PB_10x)

# Switch default assay to "RNA" for plotting of DotPlots
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results #NB: This will remove genes from the object
for (i in names (all_list)) {
  counts <- GetAssayData(all_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list[[i]] <- subset(all_list[[i]], features = rownames(counts))
}


# 4. Fix DE result lists to correct numbers of clusters and annotate them --------------------------
# 4.A DE testing - FindAllMarkers

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create a list and move DE results for correct number of clusters
FindAllMarkers_All_list <- list("all_10x" = FindAllMarkers_All_list_6$all_10x,
                                "BM_10x" = FindAllMarkers_All_list_6$BM_10x,
                                "PB_10x" = FindAllMarkers_All_list_5$PB_10x)

# Rename the clusters in the list to friendly names
FindAllMarkers_All_list$all_10x = FindAllMarkers_All_list$all_10x %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = excel_cluster_names_all_10x))
FindAllMarkers_All_list$BM_10x = FindAllMarkers_All_list$BM_10x %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = exceL_cluster_names_BM_10x))
FindAllMarkers_All_list$PB_10x = FindAllMarkers_All_list$PB_10x %>% 
  mutate(cluster = factor(cluster, levels = c(0:4), labels = excel_cluster_names_PB_10x))

# Save out the DE results
write.xlsx(FindAllMarkers_All_list,
           "Output/DE/10x/Final/FindAllMarkers_All_clusters.xlsx")

# Clean up the environment
rm(FindAllMarkers_All_list_6, FindAllMarkers_All_list_5)

# 4.A.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 

# Create a list and move DE results
FindAllMarkers_Tissue_list <- FindAllMarkers_Tissue_list_6

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list,
           "Output/DE/10x/Final/FindAllMarkers_Tissue.xlsx")

# Output Volcano plots for Tissue restricted differences
EnhancedVolcano(FindAllMarkers_Tissue_list, 
                lab = FindAllMarkers_Tissue_list$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "BM vs. PB",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 16, pointSize = 2, labSize = 8, legendLabSize = 16,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/10x/Final/Tissue/",
                "BM vs. PB", ".png"))

# Clean up the environment
rm(FindAllMarkers_Tissue_list_6, FindAllMarkers_Tissue_list_5)

# 4.A.iii DE testing - FindAllMarkers_clone - Assess global differences between clonal bins --------

# Create a list and move DE results
FindAllMarkers_Clone_list <- list("all_10x" = FindAllMarkers_Clone_list_6$all_10x,
                                  "BM_10x" = FindAllMarkers_Clone_list_6$BM_10x,
                                  "PB_10x" = FindAllMarkers_Clone_list_5$PB_10x)

# Save out the DE results
write.xlsx(FindAllMarkers_Clone_list, 
           "Output/DE/10x/Final/FindAllMarkers_Clone.xlsx")

# Clean up the environment
rm(FindAllMarkers_Clone_list_6, FindAllMarkers_Clone_list_5)

# 4.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Create a list and move DE results
FindAllMarkers_Patient_list <- list("all_10x" = FindAllMarkers_Patient_list_6$all_10x,
                                    "BM_10x" = FindAllMarkers_Patient_list_6$BM_10x,
                                    "PB_10x" = FindAllMarkers_Patient_list_5$PB_10x)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list, 
           "Output/DE/10x/Final/FindAllMarkers_Patient.xlsx")

# Clean up the environment
rm(FindAllMarkers_Patient_list_6, FindAllMarkers_Patient_list_5)

# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers Itemise differences between Clusters -------------------------------------------

# Create lists and move DE results
FindMarkers_All_all_10x_list <- FindMarkers_All_list_6[1:6]
FindMarkers_All_BM_10x_list <- FindMarkers_All_list_6[7:12]
FindMarkers_All_PB_10x_list <- FindMarkers_All_list_5[11:15]

# Rename list items to friendly names
names(FindMarkers_All_all_10x_list) <- excel_cluster_names_all_10x
names(FindMarkers_All_BM_10x_list) <- exceL_cluster_names_BM_10x
names(FindMarkers_All_PB_10x_list) <- excel_cluster_names_PB_10x

# Move all lists into another list to allow iteration for graphing
FindMarkers_All_list <- 
  list(FindMarkers_All_all_10x_list, 
       FindMarkers_All_BM_10x_list, 
       FindMarkers_All_PB_10x_list) %>%
  flatten()
names(FindMarkers_All_list) <- excel_cluster_names_full

# Save out the DE results
write.xlsx(FindMarkers_All_all_10x_list, "Output/DE/10x/Final/FindMarkers_All_all_10x.xlsx")
write.xlsx(FindMarkers_All_BM_10x_list, "Output/DE/10x/Final/FindMarkers_All_BM_10x.xlsx")
write.xlsx(FindMarkers_All_PB_10x_list, "Output/DE/10x/Final/FindMarkers_All_PB_10x.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_All_list)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_All_list[[`i`]], 
                     lab = FindMarkers_All_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_All_list)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf)) + 
    ggtitle(parse(text = cluster_names_full)[i])
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x/Final/All/",
                                   names(FindMarkers_All_list)[i], ".png"))
}

# Clean up the environment
rm(FindMarkers_All_list_6, FindMarkers_All_list_5, FindMarkers_All_all_10x_list,
   FindMarkers_All_BM_10x_list, FindMarkers_All_PB_10x_list)

# 4.B.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create a list and move DE results
FindMarkers_Tissue_list <- FindMarkers_Tissue_list_6

# Rename the items in the list
names(FindMarkers_Tissue_list) <- excel_cluster_names_all_10x

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list, 
           "Output/DE/10x/Final/FindMarkers_Tissue.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Tissue_list)) {
  temp = 
    (EnhancedVolcano(FindMarkers_Tissue_list[[`i`]], 
                     lab = FindMarkers_Tissue_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = cluster_names_all_10x[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf)) + 
    ggtitle(parse(text = cluster_names_all_10x)[i])
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x/Final/Tissue/Tissue_",
                                   names(FindMarkers_Tissue_list)[i], ".png"))
}

# Clean up the environment
rm(FindMarkers_Tissue_list_6, FindMarkers_Tissue_list_5)


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create a list and move DE results
FindConservedMarkers_Tissue_list <- FindConservedMarkers_Tissue_list_6

# Rename the items in the list
names(FindConservedMarkers_Tissue_list) <- excel_cluster_names_all_10x

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list, 
           "Output/DE/10x/Final/FindConservedMarkers_Tissue.xlsx")

# Clean up the environment
rm(FindConservedMarkers_Tissue_list_6, FindConservedMarkers_Tissue_list_5)


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAP of cluster distribution of the all_10x object
temp <- DimPlot(all_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_all_10x, 
                      labels = cluster_names_all_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "Combined") + labs(color = "Cluster") +
  scale_color_manual(values = all_10x_plot_cols, labels = label_parse()) + theme_bw() +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# Output just the cluster names for overlay purposes
temp <- DimPlot(all_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_all_10x, 
                      labels = cluster_names_all_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "Combined") + labs(color = "Cluster") +
  scale_color_manual(values = rep(c("white"),6), labels = label_parse()) +
theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
      legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "10x_clusters", ".png"), units = "cm", height = 20, width = 25)

# UMAP of cluster distribution of the BM_10x object
temp <- DimPlot(all_10x_BM, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_all_10x, 
                      labels = cluster_names_all_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "Combined - BM") + labs(color = "Cluster") +
  scale_color_manual(values = all_10x_plot_cols, labels = label_parse()) + theme_bw() + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x_BM", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the PB_10x object
temp <- DimPlot(all_10x_PB, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_all_10x, 
                      labels = cluster_names_all_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "Combined - PB") + labs(color = "Cluster") +
  scale_color_manual(values = all_10x_plot_cols, labels = label_parse()) + theme_bw() + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x_PB", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the all_10x object, grouped/split by tissue ----------------------
# Change Identity for all_10x_BM and all_10x_PB to "Tissue"
Idents(all_10x_BM) = "Tissue"
Idents(all_10x_PB) = "Tissue"

DimPlot(all_10x, group.by = "Tissue", cols = c("blue", "red")) + ggtitle(label = "Combined") + 
  labs(color = "Tissue") + theme_bw() + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x_tissue_grouped", ".png"), 
       units = "cm", height = 20, width = 25, bg = "transparent")

DimPlot(all_10x_BM, split.by = "Tissue", cols = c("blue")) + 
  ggtitle(label = "Combined - BM") + labs(color = "Tissue") + theme_bw() + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x_tissue_BM", ".png"), 
       units = "cm", height = 20, width = 25, bg = "transparent")

DimPlot(all_10x_PB, group.by = "Tissue", cols = c("red")) + 
  ggtitle(label = "Combined - PB") + labs(color = "Tissue") + theme_bw() + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "all_10x_tissue_PB", ".png"), 
       units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the BM_10x object ------------------------------------------------
temp <- DimPlot(BM_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_BM_10x, 
                      labels = cluster_names_BM_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + labs(title = "BM", color = "Cluster") +
  scale_color_manual(values = BM_10x_plot_cols, labels = label_parse()) + theme_bw() + 
  theme(legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "BM_10x", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the PB_10x object ------------------------------------------------
temp <- DimPlot(PB_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_PB_10x, 
                      labels = cluster_names_PB_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + labs(title = "PB", color = "Cluster") +
  scale_color_manual(values = PB_10x_plot_cols, labels = label_parse()) + theme_bw() + 
  theme(legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x/Final/Dimplot_",
                "PB_10x", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Gather data to plot
temp_all_10x <- as.data.frame(table(all_10x$seurat_clusters, all_10x$orig.ident))
temp_BM_10x <- as.data.frame(table(BM_10x$seurat_clusters, BM_10x$orig.ident))
temp_PB_10x <- as.data.frame(table(PB_10x$seurat_clusters, PB_10x$orig.ident))

# Convert cell numbers to proportions based on disease
temp_all_10x <- temp_all_10x %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)
temp_BM_10x <- temp_BM_10x %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)
temp_PB_10x <- temp_PB_10x %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Barplots of cluster distribution of the all_10x object -------------------------------------------
ggplot(temp_all_10x, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_all_10x), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_all_10x))) +
  labs(x = "Sample", y = "Proportion (%)", title = "Combined cluster distribution") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x/Final/Barplot_",
                "all_10x", ".png"), units = "cm", height = 20, width = 18, bg = "transparent")

# Barplots of cluster distribution of the BM_10x object --------------------------------------------
ggplot(temp_BM_10x, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_BM_10x), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_BM_10x))) +
  labs(x = "Sample", y = "Proportion (%)", title = "BM cluster distribution") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x/Final/Barplot_",
                "BM_10x", ".png"))

# Barplots of cluster distribution of the PB_10x object --------------------------------------------
ggplot(temp_PB_10x, 
       aes(fill = factor(Var1, levels = 4:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_PB_10x), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_PB_10x))) +
  labs(x = "Sample", y = "Proportion (%)", title = "PB cluster distribution") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x/Final/Barplot_",
                "PB_10x", ".png"))

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Switch clusters back to "seurat_clusters" for plotting of DotPlots - Else, dots are grey!
all_list <- map(all_list, `Idents<-`, value = "seurat_clusters")

# Extract the top 10 unique markers for the all_10x object -----------------------------------------
top_10_markers_all_10x <- 
  FindAllMarkers_All_list$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)

# Create Dotplot separated by tissue
temp = DotPlot(all_list$all_10x, 
               features = top_10_unique_all_10x, split.by = "Tissue", cols = c("blue", "red")) 

# Fix names for plot
levels(temp$data$id) = c("T[EM] ~ BM", "T[EM] ~ PB", "T[TE] ~ BM", "T[TE] ~ PB", 
                         "T[N] ~ BM", "T[N] ~ PB", "Cyto ~ T[EM] ~ BM", "Cyto ~ T[EM] ~ PB", 
                         "P[RE] ~ EX ~ BM", "P[RE] ~ EX ~ PB", "T[CM] ~ BM", "T[CM] ~ PB")

# Plot and save results
temp + theme(axis.text.x = element_text(angle = 45)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Tissue") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(hjust = 1, vjust = 1))
ggsave(filename = paste0("Output/Figures/Dotplots/10x/Final/Dotplot_",
                          "BM vs. PB", ".png"), units = "cm", height = 12, width = 34, bg = "transparent")

# Extract the top 10 unique markers for the BM_10x object ------------------------------------------
top_10_markers_BM_10x <- 
  FindAllMarkers_All_list$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM_10x <- unique(top_10_markers_BM_10x$gene)

# Create Dotplot separated by patient
temp <- DotPlot(all_list$BM_10x, 
                features = top_10_unique_BM_10x, split.by = "Patient", cols = c("green", "purple"))

# Fix names for plot
levels(temp$data$id) <- c("T[N] ~ PT43", "T[N] ~ PT63", "T[EM] ~ PT43", "T[EM] ~ PT63", 
                          "T[TE] ~ PT43", "T[TE] ~ PT63", "T[CM] ~ PT43", "T[CM] ~ PT63", 
                          "P[RE] ~ EX ~ PT43", "P[RE] ~ EX ~ PT63", 
                          "Cyto ~ T[EM] ~ PT43", "Cyto ~ T[EM] ~ PT63")

# Plot and save results
temp + theme(axis.text.x = element_text(angle = 90)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Donor - BM") + NoLegend()
ggsave(filename = paste0("Output/Figures/Dotplots/10x/Final/Dotplot_",
                          "BM pt43 vs. pt63", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB_10x object ------------------------------------------
top_10_markers_PB_10x <- 
  FindAllMarkers_All_list$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB_10x <- unique(top_10_markers_PB_10x$gene)

# Create Dotplot separated by patient
temp <- DotPlot(all_list$PB_10x, 
                features = top_10_unique_PB_10x, split.by = "Patient", 
                cols = c("#E08B00", "#0C00B8", "green", "purple"))

# Fix names for plot
levels(temp$data$id) <- c("T[TE] ~ PT13", "T[TE] ~ PT31", "T[TE] ~ PT43", "T[TE] ~ PT63", 
                          "T[EM] ~ PT13", "T[EM] ~ PT31", "T[EM] ~ PT43", "T[EM] ~ PT63", 
                          "T[CM] ~ PT13", "T[CM] ~ PT31", "T[CM] ~ PT43", "T[CM] ~ PT63", 
                          "T[N] ~ PT13", "T[N] ~ PT31", "T[N] ~ PT43", "T[N] ~ PT63", 
                          "Cyto ~ T[EM] ~ PT13", "Cyto ~ T[EM] ~ PT31", 
                          "Cyto ~ T[EM] ~ PT43", "Cyto ~ T[EM] ~ PT63")

temp + 
  theme(axis.text.x = element_text(angle = 90)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Donor - PB") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x/Final/Dotplot_",
                "PB all patients", ".png"), scale = 2:1)

# 4.D.iv Heatmap of genes that define clusters based on level of expression ------------------------

# Set Identity back to "seurat_clusters" to allow proper ordering
Idents(all_10x) = "seurat_clusters"

# Plot the mutltibarheatmap
DoMultiBarHeatmap(all_10x, features = top_10_unique_all_10x, 
                  group.by = "seurat_clusters", additional.group.by = "Tissue", label = FALSE) + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 24)) + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Heatmaps/10x/Final/Heatmap_",
                "all_10x", ".png"), scale = 2:1)


# 6. Save out data and print stats -----------------------------------------------------------------

# Move all objects into a list so we include all genes and reset default assay to integrated
all_list <- list("all_10x" = all_10x, "BM_10x" = BM_10x, "PB_10x" = PB_10x)

# Remove unnecessary cluster metadata to avoid clutter downstream
all_list$all_10x@meta.data <- all_list$all_10x@meta.data[, -c(19:42) ]
all_list$BM_10x@meta.data <- all_list$BM_10x@meta.data[, -c(19:42) ]
all_list$PB_10x@meta.data <- all_list$PB_10x@meta.data[, -c(19:42) ]

# Reinstate the integrated_snn resolution used for each object
all_list$all_10x$integrated_snn_res.0.18 <- all_10x$seurat_clusters
all_list$BM_10x$integrated_snn_res.0.28 <- BM_10x$seurat_clusters
all_list$PB_10x$integrated_snn_res.0.16 <- PB_10x$seurat_clusters

# Save the annotated Seurat objects for further analysis
saveRDS(all_list, "Data/R_out/10x/all_list.rds")

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean the environment and save out the DE analysis
rm(list = ls()[!ls() %in% c("FindAllMarkers_All_list", "FindAllMarkers_Clone_list",
                            "FindAllMarkers_Patient_list", "FindAllMarkers_Tissue_list",
                            "FindConservedMarkers_Tissue_list", "FindMarkers_All_list",
                            "FindMarkers_Tissue_list")])

# Save out the DE_results for further analysis
save.image("Data/R_out/10x/DE_Analysis_10x.rds")

gc()

print("#># Finished running '2. 10x dataset - DE Analysis' script")