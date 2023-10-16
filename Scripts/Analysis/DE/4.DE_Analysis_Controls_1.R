# Script information -------------------------------------------------------------------------------

# Title: Controls dataset #1 - DE Analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will import the previous analysis, annotate the clusters we have
# identified and create some additional figures. We will then export the Seurat objects for further
# downstream analysis. Analysis has revealed that the data is perhaps too sparse to analyse each
# condition separately and as such only the all.integrated object will be used downstream.
# Analysis has revealed the 'correct' number of clusters is 6 for the all.integrated object. We 
# will standardise the colours for clusters and ensure that plots output in a sensible manner. 
# NB: Some FindAllMarkers data will be exported as final but only the all.integrated object has been
# analysed in-depth.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '4. Controls dataset #1 - DE Analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, openxlsx, multtest, metap, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data
integrated_list <- readRDS("Data/R_out/Controls_1/integrated_list.rds")

# Load previous analysis
load("Data/R_out/Controls_1/DE_Workup_Controls_1.rds")

# Define specify_decimal function
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Define cluster names for naming of lists
cluster_names_Controls_1 <- c("TTE", "TEM", "TN", "T_eff_act", "T_eff", "TCM")

# Define character vector of plotting order of conditions
condition_order = c("all", "BM_act", "BM_rest", "PB_act", "PB_rest")

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")
# These need to be finalised
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
ccolours <- c("#B79F00", "#F8766D", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")


# 3. Pre-processing --------------------------------------------------------------------------------

# Extract the all.integrated object to allow reuse of code and so that it can be exported
all.integrated <- integrated_list$all.integrated

# Rename the clusters in the all.integrated object
all.integrated$Cluster <- all.integrated$seurat_clusters
Idents(all.integrated) <- "Cluster"
names(cluster_names_Controls_1) <- levels(all.integrated)
all.integrated <- RenameIdents(all.integrated, cluster_names_Controls_1)

# Move all objects into a list
all_list <- list("all.integrated" = all.integrated)

# Switch default assay to "RNA" for plotting of DotPlots
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results. #NB: This will remove genes from the object
for (i in names (all_list)) {
  counts <- GetAssayData(all_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list[[i]] <- subset(all_list[[i]], features = rownames(counts))
}

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create a list and move DE results
FindAllMarkers_All_list <- list("all.integrated" = FindAllMarkers_All_list_6$all.integrated)

# Rename the clusters in the list to friendly names
FindAllMarkers_All_list$all.integrated = FindAllMarkers_All_list$all.integrated %>% 
  mutate(cluster = factor(cluster, levels = c(0:5), labels = cluster_names_Controls_1))

# Save out the DE results
write.xlsx(FindAllMarkers_All_list,
           "Output/DE/Controls_1/Final/FindAllMarkers_All_clusters.xlsx")

# Clean up the environment
rm(FindAllMarkers_All_list_6, FindAllMarkers_All_list_5)

# 4.A.ii DE testing - FindAllMarkers_Stim - Assess global differences between Stimulation state ----

# Create a list and move DE results
FindAllMarkers_Stim_list <- FindAllMarkers_Stim_list_6

# Save out the DE results
write.xlsx(FindAllMarkers_Stim_list,
           "Output/DE/Controls_1/Final/FindAllMarkers_Stim.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindAllMarkers_Stim_list)) {
  temp <- 
    (EnhancedVolcano(FindAllMarkers_Stim_list[[`i`]], 
                     lab = FindAllMarkers_Stim_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindAllMarkers_Stim_list)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Final/FindAllMarkers_Stim/",
                                   names(FindAllMarkers_Stim_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindAllMarkers_Stim_list_6, FindAllMarkers_Stim_list_5)

# 4.A.iii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------

# Create a list and move DE results
FindAllMarkers_Tissue_list <- FindAllMarkers_Tissue_list_6

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list,
           "Output/DE/Controls_1/Final/FindAllMarkers_Tissue.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindAllMarkers_Tissue_list)) {
  temp <- 
    (EnhancedVolcano(FindAllMarkers_Tissue_list[[`i`]], 
                     lab = FindAllMarkers_Tissue_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindAllMarkers_Tissue_list)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Final/FindAllMarkers_Tissue/",
                                   names(FindAllMarkers_Tissue_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindAllMarkers_Tissue_list_6, FindAllMarkers_Tissue_list_5)

# 4.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Create a list and move DE results
FindAllMarkers_Patient_list <- FindAllMarkers_Patient_list_6$all.integrated

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list, 
           "Output/DE/Controls_1/Final/FindAllMarkers_Patient.xlsx")

# Clean up the environment
rm(FindAllMarkers_Patient_list_6, FindAllMarkers_Patient_list_5)


# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers Itemise differences between Clusters -------------------------------------------

# Create lists and move DE results
FindMarkers_All_list <- FindMarkers_All_list_6[25:30]
# Rename list items to friendly names
names(FindMarkers_All_list) <- cluster_names_Controls_1

# Save out the DE results
write.xlsx(FindMarkers_All_list, 
           "Output/DE/Controls_1/Final/FindMarkers_All.xlsx")

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
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Final/Cluster/",
                                   names(FindMarkers_All_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_All_list_6, FindMarkers_All_list_5)

# 4.B.ii FindMarkers - Itemise differences between activated and resting ---------------------------

# Create a list and move DE results
FindMarkers_Stim_list <- FindMarkers_Stim_list_6[13:18]

# Rename the items in the list
names(FindMarkers_Stim_list) <- cluster_names_Controls_1

# Save out the DE results
write.xlsx(FindMarkers_Stim_list, 
           "Output/DE/Controls_1/Final/FindMarkers_Stim.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Stim_list)) {
  temp = 
    (EnhancedVolcano(FindMarkers_Stim_list[[`i`]], 
                     lab = FindMarkers_Stim_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Stim_list)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 6, 
                     legendPosition = "none", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Final/FindMarkers_Stim/",
                                   names(FindMarkers_Stim_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_Stim_list_6, FindMarkers_Stim_list_5)

# 4.B.iii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create a list and move DE results
FindMarkers_Tissue_list <- FindMarkers_Tissue_list_6[7:12]

# Rename the items in the list
names(FindMarkers_Tissue_list) <- cluster_names_Controls_1

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list, 
           "Output/DE/Controls_1/Final/FindMarkers_Tissue.xlsx")

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
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Final/FindMarkers_Tissue/",
                                   names(FindMarkers_Tissue_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_Tissue_list_6, FindMarkers_Tissue_list_5)


# 4.C FindConservedMarkers - Find conserved markers across condition - Itemised --------------------

# 4.C.i FindConservedMarkers - Averages across Stimulation status ----------------------------------

# Create a list and move DE results
FindConservedMarkers_Stim_list <- FindConservedMarkers_Stim_list_6[13:18]

# Rename the items in the list
names(FindConservedMarkers_Stim_list) <- cluster_names_Controls_1

# Save out the DE results
write.xlsx(FindConservedMarkers_Stim_list, 
           "Output/DE/Controls_1/Final/FindConservedMarkers_Stim.xlsx")

# Clean up the environment
rm(FindConservedMarkers_Stim_list_6, FindConservedMarkers_Stim_list_5)

# 4.C.ii FindConservedMarkers - Averages across Tissue ---------------------------------------------

# Create a list and move DE results
FindConservedMarkers_Tissue_list <- FindConservedMarkers_Tissue_list_6[7:12]

# Rename the items in the list
names(FindConservedMarkers_Tissue_list) <- cluster_names_Controls_1

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list, 
           "Output/DE/Controls_1/Final/FindConservedMarkers_Tissue.xlsx")

# Clean up the environment
rm(FindConservedMarkers_Tissue_list_6, FindConservedMarkers_Tissue_list_5)


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAP of cluster distribution of the all.integrated object
DimPlot(all.integrated, label = TRUE, cols = ccolours) + 
  labs(title = "Controls_1", color = "Cluster") +
  theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Dimplots/Controls_1/Final/Dimplot_",
                "Controls_1", ".png"))

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Create barplot of distribution by stimulation status
# Gather data to plot
temp_Stim <- 
  as.data.frame(table(all.integrated$seurat_clusters, all.integrated$stimulation_status))
# Convert cell numbers to proportions
temp_Stim <- temp_Stim %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells
ggplot(temp_Stim, 
       aes(fill = factor(Var1, levels = 5:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  labs(x = "Stimulation status", y = "Cells", 
       title = "Controls_1 - Cluster distribution by stimulation status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_cluster_by_stimulation_status", ".png"))

# Plot as proportion
ggplot(temp_Stim, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  labs(x = "Stimulation status", y = "Proportion (%)", 
       title = "Controls_1 - Cluster distribution by stimulation status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_cluster_by_stimulation_status_proportions", ".png"))

# Create barplot of distribution by Tissue
# Gather data to plot
temp_Tissue <- 
  as.data.frame(table(all.integrated$seurat_clusters, all.integrated$Tissue))
# Convert cell numbers to proportions
temp_Tissue <- temp_Tissue %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells
ggplot(temp_Tissue, 
       aes(fill = factor(Var1, levels = 5:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  labs(x = "Tissue", y = "Cells", 
       title = "Controls_1 - Cluster distribution by Tissue") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_cluster_by_Tissue", ".png"))

# Plot as proportion
ggplot(temp_Tissue, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  labs(x = "Tissue", y = "Proportion (%)", 
       title = "Controls_1 - Cluster distribution by Tissue") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_cluster_by_Tissue_proportions", ".png"))

# Calculate distribution of clusters on a per sample level
temp_all.integrated <- 
  as.data.frame(table(all.integrated$seurat_clusters, all.integrated$orig.ident))
# Convert cell numbers to proportions
temp_all.integrated <- 
  temp_all.integrated %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells
ggplot(temp_all.integrated, 
       aes(fill = factor(Var1, levels = 5:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  geom_text(aes(label = paste(Freq)), position = position_stack(vjust = 0.5)) +
  labs(x = "Donor", y = "Cells", 
       title = "Controls_1 - Cluster distrubution by sample") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_clusters_by_sample", ".png"))

# Plot as proportion
ggplot(temp_all.integrated, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = cluster_names_Controls_1) + 
  geom_text(aes(label = paste(specify_decimal(Proportion, 2))), position = position_stack(vjust = 0.5)) +
  labs(x = "Donor", y = "Proportion", 
       title = "Controls_1 - Cluster distrubution by sample") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_clusters_by_sample_proportions", ".png"))

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all.integrated object
top_10_markers_all.integrated <- 
  FindAllMarkers_All_list$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all.integrated <- unique(top_10_markers_all.integrated$gene)

# Plot as a Dotplot separated by stimulation status
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all.integrated, 
        split.by = "stimulation_status", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Stimulation status") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Final/Dotplot_",
                "Activated vs. Resting", ".png"), scale = 2:1)

# Plot as a Dotplot separated by tissue
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all.integrated, 
        split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Tissue") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Final/Dotplot_",
                "BM vs. PB", ".png"), scale = 2:1)

# Plot as a Dotplot separated by patient
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all.integrated, split.by = "orig.ident", 
        cols = c("blue", "red")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "Gene expression by Donor") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Final/Dotplot_",
                "Donor1 vs. Donor2", ".png"), scale = 2:1)

# 4.D.iv Heatmap of genes that define clusters based on level of expression ------------------------

# Delete the barcode metadata - not sure why this is required but OK
#temp = all.integrated
#temp$barcode = NULL
#
# Plot MulitHeatMap
#png("Output/Figures/Seurat/Controls_1/Final/Heatmap/Heatmap_all.integrated.png")
#plot_heatmap(temp,
#             n = 10,
#             markers = top_10_unique_all.integrated,
#             sort_var = c("Cluster","Tissue", "stimulation_status"),
#             anno_var = c("Cluster", "Tissue", "stimulation_status"),
#             anno_colors = list(ccolours,
#                                c("blue","red"), 
#                                c("green", "purple")),
#             hm_limit = c(-2,0,2),
#             hm_colors = c("purple","black","yellow"))
#dev.off()

# 5. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Save the annotated Seurat objects for further analysis
saveRDS(all.integrated, "Data/R_out/Controls_1/Controls_1_complete.rds")

# Clean the environment and save out the DE analysis
rm(list = ls()[!ls() %in% c("FindAllMarkers_All_list", "FindAllMarkers_Patient_list",
                            "FindAllMarkers_Stim_list", "FindAllMarkers_Tissue_list",
                            "FindConservedMarkers_Stim_list", "FindConservedMarkers_Tissue_list",
                            "FindMarkers_All_list", "FindMarkers_Stim_list", 
                            "FindMarkers_Tissue_list")])

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_1/DE_Analysis_Controls_1.rds")

print("#># Finished running '4. Controls dataset #1 - DE Analysis' script")