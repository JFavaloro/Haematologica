# Script information -------------------------------------------------------------------------------

# Title: ProjecTILs/10x Reference atlases - DE Analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will import the previous analysis, annotate the clusters we have
# identified and create some additional figures. We will then export the reference atlases for
# downstream analysis. Analysis has revealed the 'correct' number of clusters is 6 for all atlases.
# We will standardise the colours for clusters and ensure that plots output in a sensible manner.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. ProjecTILs/10x Reference atlases - DE Analysis' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, openxlsx, multtest, metap, EnhancedVolcano, scales, DoMultiBarHeatmap)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load reference atlases
all_refs <- readRDS("Data/R_out/ProjecTILs/all_refs.rds")

# Load previous analysis
load("Data/R_out/ProjecTILs/DE_Workup_10x_atlases.rds")

# Create character vectors for naming of clusters and lists
cluster_names_full <- 
  c("Combined ~ T[EM]", "Combined ~ T[TE]", "Combined ~ T[N]", 
    "Combined ~ IL7R^`+` ~ T[M]", "Combined ~ Cyto ~ T[EM]", "Combined ~ Activated",
    "BM ~ T[EM]", "BM ~ T[N]", "BM ~ IL7R^`+` ~ T[M]", 
    "BM ~ T[TE]", "BM ~ Activated", "BM ~ Cyto ~ T[EM]", 
    "PB ~ T[TE]", "PB ~ T[EM]", "PB ~ IL7R^`+` ~ T[M]", 
    "PB ~ T[N]", "PB ~ Activated", "PB ~ IFN")

excel_cluster_names_full <- 
  c("all_10x_atlas_TEM", "all_10x_atlas_TTE", "all_10x_atlas_TN", "all_10x_atlas_IL7R_TM", 
    "all_10x_atlas_Cyto_TEM", "all_10x_atlas_Activated",
    "BM_10x_atlas_TEM", "BM_10x_atlas_TN", "BM_10x_atlas_IL7R_TM", "BM_10x_atlas_TTE", 
    "BM_10x_atlas_Activated", "BM_10x_atlas_Cyto_TEM", 
    "PB_10x_atlas_TTE", "PB_10x_atlas_TEM", "PB_10x_atlas_IL7R_TM", "PB_10x_atlas_TN", 
    "PB_10x_atlas_Activated", "PB_10x_atlas_IFN")

cluster_names_all_10x <- c("T[EM]", "T[TE]", "T[N]", "IL7R^`+` ~ T[M]", "Cyto ~ T[EM]", "Activated")
cluster_names_BM_10x <- c("T[EM]", "T[N]", "IL7R^`+` ~ T[M]", "T[TE]", "Activated", "Cyto ~ T[EM]")
cluster_names_PB_10x <- c("T[TE]", "T[EM]", "IL7R^`+` ~ T[M]", "T[N]", "Activated", "IFN")

excel_cluster_names_all_10x <- c("TEM", "TTE", "TN", "IL7R_TM", "Cyto_TEM", "Activated")
excel_cluster_names_BM_10x <- c("TEM", "TN", "IL7R_TM", "TTE", "Activated", "Cyto_TEM")
excel_cluster_names_PB_10x <- c("TTE", "TEM", "IL7R_TM", "TN", "Activated", "IFN")

# Create character vectors of colours for plotting clusters
all_10x_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[TE]" = "#B79F00", "T[N]" = "#00BA38", 
    "IL7R^`+` ~ T[M]" = "#F564E3", "Cyto ~ T[EM]" = "#00BFC4", "Activated" = "#8A2BE2")

BM_10x_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[N]" = "#00BA38", "IL7R^`+` ~ T[M]" = "#F564E3", 
    "T[TE]" = "#B79F00", "Activated" = "#8A2BE2", "Cyto ~ T[EM]" = "#00BFC4")

PB_10x_plot_cols <- 
  c("T[TE]" = "#B79F00", "T[EM]" = "#F8766D", "IL7R^`+` ~ T[M]" = "#F564E3", 
    "T[N]" = "#00BA38", "Activated" = "#8A2BE2", "IFN" = "#FF69B4")

ccolours_all_10x <- c("#F8766D", "#B79F00", "#00BA38", "#F564E3", "#00BFC4", "#8A2BE2")
ccolours_BM_10x <- c("#F8766D", "#00BA38", "#F564E3", "#B79F00", "#8A2BE2", "#00BFC4")
ccolours_PB_10x <- c("#B79F00", "#F8766D", "#F564E3", "#00BA38", "#8A2BE2", "#FF69B4")


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Extract the reference objects we're interested in and set the correct number of clusters -----
all_10x_atlas <- all_refs$all_10x
BM_10x_atlas <- all_refs$BM_10x
PB_10x_atlas <- all_refs$PB_10x

# 3.2 Add Sample metadata --------------------------------------------------------------------------
all_10x_atlas[["Sample"]] = "MM Reference (BM+PB)"
BM_10x_atlas[["Sample"]] = "MM Reference (BM)"
PB_10x_atlas[["Sample"]] = "MM Reference (PB)"


# 3.3 Set resolutions of Seurat objects to allow clustering into 6 clusters ------------------------
# Set resolution to 0.32 for all_10x_atlas object
all_10x_atlas$seurat_clusters <- 
  all_10x_atlas$integrated_snn_res.0.32
Idents(all_10x_atlas) <- "seurat_clusters"

# Set resolution to 0.36 for BM_10x_atlas object
BM_10x_atlas$seurat_clusters <- 
  BM_10x_atlas$integrated_snn_res.0.36
Idents(BM_10x_atlas) <- "seurat_clusters"

# Set resolution to 0.25 for PB_10x_atlas object
PB_10x_atlas$seurat_clusters <- 
  PB_10x_atlas$integrated_snn_res.0.25
Idents(PB_10x_atlas) <- "seurat_clusters"

# 3.4.i Rename the all_10x atlas object and set the atlas Idents accordingly -----------------------
f.clusters <- rep(NA, dim(all_10x_atlas)[2])
names(f.clusters) <- colnames(all_10x_atlas)
f.clusters[all_10x_atlas$seurat_clusters %in% c(0)] = "T[EM]"
f.clusters[all_10x_atlas$seurat_clusters %in% c(1)] = "T[TE]"
f.clusters[all_10x_atlas$seurat_clusters %in% c(2)] = "T[N]"
f.clusters[all_10x_atlas$seurat_clusters %in% c(3)] = "IL7R^`+` ~ T[M]"
f.clusters[all_10x_atlas$seurat_clusters %in% c(4)] = "Cyto ~ T[EM]"
f.clusters[all_10x_atlas$seurat_clusters %in% c(5)] = "Activated"
all_10x_atlas <- AddMetaData(all_10x_atlas, as.factor(f.clusters), col.name = "functional.cluster")
Idents(all_10x_atlas) <- "functional.cluster"
levels(all_10x_atlas) <- cluster_names_all_10x
all_10x_atlas$functional.cluster <- 
  fct_relevel(all_10x_atlas$functional.cluster, cluster_names_all_10x)

# 3.4.ii Rename the BM 10x atlas object and set the atlas Idents accordingly -----------------------
f.clusters <- rep(NA, dim(BM_10x_atlas)[2])
names(f.clusters) <- colnames(BM_10x_atlas)
f.clusters[BM_10x_atlas$seurat_clusters %in% c(0)] = "T[EM]"
f.clusters[BM_10x_atlas$seurat_clusters %in% c(1)] = "T[N]"
f.clusters[BM_10x_atlas$seurat_clusters %in% c(2)] = "IL7R^`+` ~ T[M]"
f.clusters[BM_10x_atlas$seurat_clusters %in% c(3)] = "T[TE]"
f.clusters[BM_10x_atlas$seurat_clusters %in% c(4)] = "Activated"
f.clusters[BM_10x_atlas$seurat_clusters %in% c(5)] = "Cyto ~ T[EM]"
BM_10x_atlas <- AddMetaData(BM_10x_atlas, as.factor(f.clusters), col.name = "functional.cluster")
Idents(BM_10x_atlas) <- "functional.cluster"
levels(BM_10x_atlas) <- cluster_names_BM_10x
BM_10x_atlas$functional.cluster <- 
  fct_relevel(BM_10x_atlas$functional.cluster, cluster_names_BM_10x)

# 3.4.iii Rename the PB 10x atlas object and set the atlas Idents accordingly ----------------------
f.clusters <- rep(NA, dim(PB_10x_atlas)[2])
names(f.clusters) <- colnames(PB_10x_atlas)
f.clusters[PB_10x_atlas$seurat_clusters %in% c(0)] = "T[TE]"
f.clusters[PB_10x_atlas$seurat_clusters %in% c(1)] = "T[EM]"
f.clusters[PB_10x_atlas$seurat_clusters %in% c(2)] = "IL7R^`+` ~ T[M]"
f.clusters[PB_10x_atlas$seurat_clusters %in% c(3)] = "T[N]"
f.clusters[PB_10x_atlas$seurat_clusters %in% c(4)] = "Activated"
f.clusters[PB_10x_atlas$seurat_clusters %in% c(5)] = "IFN"
PB_10x_atlas <- AddMetaData(PB_10x_atlas, as.factor(f.clusters), col.name = "functional.cluster")
Idents(PB_10x_atlas) <- "functional.cluster"
levels(PB_10x_atlas) <- cluster_names_PB_10x
PB_10x_atlas$functional.cluster <- 
  fct_relevel(PB_10x_atlas$functional.cluster, cluster_names_PB_10x)


# 3.4 Continue Pre-processing ----------------------------------------------------------------------

# Move all objects into a list
all_list <- list("all_10x" = all_10x_atlas, "BM_10x" = BM_10x_atlas, "PB_10x" = PB_10x_atlas)

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
FindAllMarkers_All_list <- list("all_10x_atlas" = FindAllMarkers_All_list_6$all_10x,
                                "BM_10x_atlas" = FindAllMarkers_All_list_6$BM_10x,
                                "PB_10x_atlas" = FindAllMarkers_All_list_6$PB_10x)

# Rename the clusters in the list to friendly names
FindAllMarkers_All_list$all_10x_atlas = FindAllMarkers_All_list$all_10x_atlas %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = excel_cluster_names_all_10x))
FindAllMarkers_All_list$BM_10x_atlas = FindAllMarkers_All_list$BM_10x_atlas %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = excel_cluster_names_BM_10x))
FindAllMarkers_All_list$PB_10x_atlas = FindAllMarkers_All_list$PB_10x_atlas %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = excel_cluster_names_PB_10x))

# Save out the DE results
write.xlsx(FindAllMarkers_All_list,
           "Output/DE/10x_atlases/Final/FindAllMarkers_All_10x_atlases.xlsx")

# Clean up the environment
rm(FindAllMarkers_All_list_7, FindAllMarkers_All_list_6)

# 4.A.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Create a list and move DE results
FindAllMarkers_Tissue_list <- FindAllMarkers_Tissue_list_6

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list,
           "Output/DE/10x_atlases/Final/FindAllMarkers_Tissue_all_10x_atlas.xlsx")

# Output Volcano plots for Tissue restricted differences
EnhancedVolcano(FindAllMarkers_Tissue_list, 
                lab = FindAllMarkers_Tissue_list$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "BM vs. PB",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/10x_atlases/Final/Tissue/",
                "BM vs. PB", ".png"), scale = 2:1)

# Clean up the environment
rm(FindAllMarkers_Tissue_list_7, FindAllMarkers_Tissue_list_6)

# 4.A.iii DE testing - FindAllMarkers_clone - Assess global differences between clonal bins --------

# Create a list and move DE results
FindAllMarkers_Clone_list <- list("all_10x_atlas" = FindAllMarkers_Clone_list_6$all_10x,
                                  "BM_10x_atlas" = FindAllMarkers_Clone_list_6$BM_10x,
                                  "PB_10x_atlas" = FindAllMarkers_Clone_list_6$PB_10x)

# Save out the DE results
write.xlsx(FindAllMarkers_Clone_list, 
           "Output/DE/10x_atlases/Final/FindAllMarkers_Clone_10x_atlases.xlsx")

# Clean up the environment
rm(FindAllMarkers_Clone_list_7, FindAllMarkers_Clone_list_6)

# 4.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Create a list and move DE results
FindAllMarkers_Patient_list <- list("all_10x_atlas" = FindAllMarkers_Patient_list_6$all_10x,
                                    "BM_10x_atlas" = FindAllMarkers_Patient_list_6$BM_10x,
                                    "PB_10x_atlas" = FindAllMarkers_Patient_list_6$PB_10x)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list, 
           "Output/DE/10x_atlases/Final/FindAllMarkers_Patient_10x_atlases.xlsx")

# Clean up the environment
rm(FindAllMarkers_Patient_list_7, FindAllMarkers_Patient_list_6)

# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers Itemise differences between Clusters -------------------------------------------

# Create lists and move DE results
FindMarkers_All_all_10x_list <- FindMarkers_All_list_6[1:6]
FindMarkers_All_BM_10x_list <- FindMarkers_All_list_6[7:12]
FindMarkers_All_PB_10x_list <- FindMarkers_All_list_6[13:18]

# Rename list items to friendly names
names(FindMarkers_All_all_10x_list) <- excel_cluster_names_all_10x
names(FindMarkers_All_BM_10x_list) <- excel_cluster_names_BM_10x
names(FindMarkers_All_PB_10x_list) <- excel_cluster_names_PB_10x

# Move all lists into another list to allow iteration for graphing
FindMarkers_All_list <- 
  list(FindMarkers_All_all_10x_list, 
       FindMarkers_All_BM_10x_list, 
       FindMarkers_All_PB_10x_list) %>%
  flatten()
names(FindMarkers_All_list) <- excel_cluster_names_full

# Save out the DE results
write.xlsx(FindMarkers_All_all_10x_list, 
           "Output/DE/10x_atlases/Final/FindMarkers_All_all_10x_atlas.xlsx")
write.xlsx(FindMarkers_All_BM_10x_list, 
           "Output/DE/10x_atlases/Final/FindMarkers_All_BM_10x_atlas.xlsx")
write.xlsx(FindMarkers_All_PB_10x_list, 
           "Output/DE/10x_atlases/Final/FindMarkers_All_PB_10x_atlas.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_All_list)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_All_list[[`i`]], 
                     lab = FindMarkers_All_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_All_list)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf)) + 
    ggtitle(parse(text = cluster_names_full)[i])
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Final/All/",
                                   names(FindMarkers_All_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_All_list_7, FindMarkers_All_list_6, FindMarkers_All_all_10x_list,
   FindMarkers_All_BM_10x_list, FindMarkers_All_PB_10x_list)

# 4.B.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create a list and move DE results
FindMarkers_Tissue_list <- FindMarkers_Tissue_list_6

# Rename the items in the list
names(FindMarkers_Tissue_list) <- excel_cluster_names_all_10x

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list, 
           "Output/DE/10x_atlases/Final/FindMarkers_Tissue_all_10x_atlas.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Tissue_list)) {
  temp = 
    (EnhancedVolcano(FindMarkers_Tissue_list[[`i`]], 
                     lab = FindMarkers_Tissue_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Tissue_list)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf)) + 
    ggtitle(parse(text = cluster_names_all_10x)[i])
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Final/Tissue/",
                                   names(FindMarkers_Tissue_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_Tissue_list_7, FindMarkers_Tissue_list_6)


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create a list and move DE results
FindConservedMarkers_Tissue_list <- FindConservedMarkers_Tissue_list_6

# Rename the items in the list
names(FindConservedMarkers_Tissue_list) <- excel_cluster_names_all_10x

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list, 
           "Output/DE/10x_atlases/Final/FindConservedMarkers_Tissue_all_10x_atlas.xlsx")

# Clean up the environment
rm(FindConservedMarkers_Tissue_list_7, FindConservedMarkers_Tissue_list_6)

# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAP of cluster distribution of the all_10x object
temp <- DimPlot(all_10x_atlas, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_all_10x, 
                      labels = cluster_names_all_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "Combined") + labs(color = "Cluster") +
  scale_color_manual(values = all_10x_plot_cols, labels = label_parse()) + theme_bw() +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x_atlases/Final/Dimplot_",
                "all_10x_atlas", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the BM_10x object ------------------------------------------------
temp <- DimPlot(BM_10x_atlas, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_BM_10x, 
                      labels = cluster_names_BM_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "MM BM reference atlas") + 
  labs(color = "Cluster") +
  scale_color_manual(values = BM_10x_plot_cols, labels = label_parse()) + theme_bw() +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x_atlases/Final/Dimplot_",
                "BM_10x_atlas", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# UMAP of cluster distribution of the PB_10x object ------------------------------------------------
temp <- DimPlot(PB_10x_atlas, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp <- LabelClusters(plot = temp, id = "ident", clusters = cluster_names_PB_10x, 
                      labels = cluster_names_PB_10x, parse = TRUE, size = 10)
temp + guides(fill = guide_legend("Cluster")) + ggtitle(label = "MM PB reference atlas") + 
  labs(color = "Cluster") +
  scale_color_manual(values = PB_10x_plot_cols, labels = label_parse()) + theme_bw() +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Dimplots/10x_atlases/Final/Dimplot_",
                "PB_10x_atlas", ".png"), units = "cm", height = 20, width = 25, bg = "transparent")

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Gather data to plot
temp_all_10x <- as.data.frame(table(all_10x_atlas$seurat_clusters, all_10x_atlas$orig.ident))
temp_BM_10x <- as.data.frame(table(BM_10x_atlas$seurat_clusters, BM_10x_atlas$orig.ident))
temp_PB_10x <- as.data.frame(table(PB_10x_atlas$seurat_clusters, PB_10x_atlas$orig.ident))

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
  labs(x = "Sample", y = "Proportion (%)", title = "Combined reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x_atlases/Final/Barplot_",
                "all_10x_atlas", ".png"), units = "cm", height = 20, width = 18, bg = "transparent")

# Barplots of cluster distribution of the BM_10x object --------------------------------------------
ggplot(temp_BM_10x, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_BM_10x), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_BM_10x))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "MM BM reference atlas cluster distribution") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x_atlases/Final/Barplot_",
                "BM_10x_atlas", ".png"), units = "cm", height = 20, width = 18, bg = "transparent")

# Barplots of cluster distribution of the PB_10x object --------------------------------------------
ggplot(temp_PB_10x, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_PB_10x), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_PB_10x))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "MM PB reference atlas cluster distribution") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/10x_atlases/Final/Barplot_",
                "PB_10x_atlas", ".png"), units = "cm", height = 20, width = 18, bg = "transparent")

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Rename Idents to allow plotting with colours
levels(all_list$all_10x$functional.cluster) <- c("A", "B", "C", "D", "E", "F")
levels(all_list$BM_10x$functional.cluster) <- c("A", "B", "C", "D", "E", "F")
levels(all_list$PB_10x$functional.cluster) <- c("A", "B", "C", "D", "E", "F")
all_list <- map(all_list, `Idents<-`, value = "functional.cluster")

# Extract the top 10 unique markers for the all_10x object -----------------------------------------
top_10_markers_all_10x <- 
  FindAllMarkers_All_list$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)

# Create Dotplot separated by tissue
temp <- DotPlot(all_list$all_10x, 
                features = top_10_unique_all_10x, split.by = "Tissue", cols = c("blue", "red"))

# Fix names for plot
levels(temp$data$id) = c("T[EM] ~ BM", "T[EM] ~ PB", "T[TE] ~ BM", "T[TE] ~ PB", 
                         "T[N] ~ BM", "T[N] ~ PB", "IL7R^`+` ~ T[M] ~ BM", "IL7R^`+` ~ T[M] ~ PB", 
                         "Cyto ~ T[EM] ~ BM", "Cyto ~ T[EM] ~ PB", 
                         "Activated ~ BM", "Activated ~ PB")

# Plot and save results
temp + theme(axis.text.x = element_text(angle = 90)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", 
       title = "Gene expression by Tissue - MM BM+PB reference atlas") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Final/Dotplot_",
                "BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM_10x object ------------------------------------------
top_10_markers_BM_10x <- 
  FindAllMarkers_All_list$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM_10x <- unique(top_10_markers_BM_10x$gene)

# Create Dotplot separated by patient
temp <- DotPlot(all_list$BM_10x, 
                features = top_10_unique_BM_10x, split.by = "Patient", cols = c("green", "purple"))

# Fix names for plot
levels(temp$data$id) <- c("T[EM] ~ PT43", "T[EM] ~ PT63", "T[N] ~ PT43", "T[N] ~ PT63", 
                          "IL7R^`+` ~ T[M] ~ PT43", "IL7R^`+` ~ T[M] ~ PT63", "T[TE] ~ PT43", "T[TE] ~ PT63", 
                          "Activated ~ PT43", "Activated ~ EX ~ PT63", 
                          "Cyto ~ T[EM] ~ PT43", "Cyto ~ T[EM] ~ PT63")

# Plot and save results
temp + theme(axis.text.x = element_text(angle = 90)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", 
       title = "Gene expression by Donor - MM BM reference atlas") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Final/Dotplot_",
                "BM pt43 vs. pt63", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB_10x object ------------------------------------------
top_10_markers_PB_10x <- 
  FindAllMarkers_All_list$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB_10x <- unique(top_10_markers_PB_10x$gene)

# Create Dotplot separated by patient
temp <- DotPlot(all_list$PB_10x, 
                features = top_10_unique_PB_10x, split.by = "Patient", 
                cols = c("orange", "turquoise", "green", "purple")) 

# Fix names for plot
levels(temp$data$id) <- c("T[TE] ~ PT13", "T[TE] ~ PT31", "T[TE] ~ PT43", "T[TE] ~ PT63", 
                          "T[EM] ~ PT13", "T[EM] ~ PT31", "T[EM] ~ PT43", "T[EM] ~ PT63", 
                          "IL7R^`+` ~ T[M] ~ PT13", "IL7R^`+` ~ T[M] ~ PT31", "IL7R^`+` ~ T[M] ~ PT43", "IL7R^`+` ~ T[M] ~ PT63", 
                          "T[N] ~ PT13", "T[N] ~ PT31", "T[N] ~ PT43", "T[N] ~ PT63", 
                          "Activated ~ PT13", "Activated ~ PT31", 
                          "Activated ~ PT43", "Activated ~ PT63",
                          "IFN ~ PT13", "IFN ~ PT31", 
                          "IFN ~ PT43", "IFN ~ PT63")

# Plot and save results
temp + theme(axis.text.x = element_text(angle = 90)) + scale_y_discrete(labels = label_parse()) +
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Donor - PB_10x_atlas") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Final/Dotplot_",
                "PB all patients", ".png"), scale = 2:1)

# 4.D.iv Heatmap of genes that define clusters based on level of expression ------------------------

temp <- DoHeatmap(BM_10x_atlas, features = top_10_unique_BM_10x, assay = "integrated", 
                  group.bar = TRUE, group.colors = ccolours_BM_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp + ggtitle(label = "MM BM reference atlas") +
  labs(color = "Cluster", 
       labels = label_parse()) + scale_color_manual(values = BM_10x_plot_cols) + 
  theme(text=element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Heatmaps/10x_atlases/Final/",
                "MM BM reference atlas", ".png"), units = "cm", height = 30, width = 50, bg = "transparent")

temp <- DoHeatmap(PB_10x_atlas, features = top_10_unique_PB_10x, assay = "integrated", 
                  group.bar = TRUE, group.colors = ccolours_PB_10x, label = FALSE, combine = FALSE)
temp <- temp[[1]]
temp + ggtitle(label = "MM PB reference atlas") +
  labs(color = "Cluster", 
       labels = label_parse()) + scale_color_manual(values = PB_10x_plot_cols) + 
  theme(text=element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.text.align = 0)
ggsave(filename = 
         paste0("Output/Figures/Heatmaps/10x_atlases/Final/",
                "MM PB reference atlas", ".png"), units = "cm", height = 30, width = 50, bg = "transparent")


# 5. Save out data and print stats -----------------------------------------------------------------

# Remove unnecessary cluster metadata to avoid clutter downstream
all_10x_atlas@meta.data <- all_10x_atlas@meta.data[, -c(27:38) ]
BM_10x_atlas@meta.data <- BM_10x_atlas@meta.data[, -c(27:37) ]
PB_10x_atlas@meta.data <- PB_10x_atlas@meta.data[, -c(27:38) ]

# Save the 10x atlases for downstream use
saveRDS(all_10x_atlas, "Data/R_out/ProjecTILs/all_10x_atlas.rds")
saveRDS(BM_10x_atlas, "Data/R_out/ProjecTILs/BM_10x_atlas.rds")
saveRDS(PB_10x_atlas, "Data/R_out/ProjecTILs/PB_10x_atlas.rds")

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean the environment and save out the DE analysis
rm(list = ls()[!ls() %in% c("FindAllMarkers_All_list", "FindAllMarkers_Clone_list",
                            "FindAllMarkers_Patient_list", "FindAllMarkers_Tissue_list",
                            "FindConservedMarkers_Tissue_list", "FindMarkers_All_list",
                            "FindMarkers_Tissue_list")])

# Save out the DE_results for further analysis
save.image("Data/R_out/ProjecTILs/DE_Analysis_10x_atlases.rds")

gc()

print("#># Finished running '5. ProjecTILs/10x Reference atlases - DE Analysis' script")