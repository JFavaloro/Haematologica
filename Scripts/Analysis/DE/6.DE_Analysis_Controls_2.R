# Script information -------------------------------------------------------------------------------

# Title: Controls dataset #2 - DE Analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will import the previous analysis, annotate the clusters we have
# identified and create some additional figures. We will then export the Seurat objects for further
# downstream analysis. Analysis has revealed the 'correct' number of clusters is 5 for this dataset.
# We will standardise the colours for clusters and ensure that plots output in a sensible manner.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '6. Controls dataset #2 - DE Analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, openxlsx, multtest, metap, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data
integrated_list <- readRDS("Data/R_out/Controls_2/integrated_list.rds")

# Load previous analysis
load("Data/R_out/Controls_2/DE_Workup_Controls_2.rds")

# Define cluster names for naming of lists
cluster_names_Controls_2 <- c("TCM", "TN", "TTE", "TEM", "Cyto_TEM")

# Create character vector of plotting order for disease
disease_order = c("HD", "MGUS", "SMM", "MM")

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")
# These need to be finalised
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
ccolours<- c("#F564E3", "#00BA38", "#B79F00", "#F8766D", "#619CFF")


# 3. Pre-processing --------------------------------------------------------------------------------

# Extract the all.integrated object to allow reuse of code and so that it can be exported
all.integrated <- integrated_list$all.integrated

# Set resolution to 0.20 for all.integrated object to allow the correct number of clusters (5)
all.integrated$seurat_clusters <- all.integrated$integrated_snn_res.0.2
Idents(all.integrated) <- "seurat_clusters"

# Rename the clusters in the all.integrated object
all.integrated$Cluster <- all.integrated$seurat_clusters
Idents(all.integrated) <- "Cluster"
names(cluster_names_Controls_2) <- levels(all.integrated)
all.integrated <- RenameIdents(all.integrated, cluster_names_Controls_2)

# Change the order of disease for more sensible results
all.integrated$Disease <- factor(all.integrated$Disease, levels = disease_order)

# Change the order of orig.ident for more sensible results
all.integrated$orig.ident <- 
  factor(all.integrated$orig.ident, 
         levels = c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7", "HD8", "HD9", 
                    "MGUS1", "MGUS2", "MGUS3", "MGUS4", "MGUS5", 
                    "SMM1", "SMM2", "SMM3", "SMM4", "SMM5", "SMM6", "SMM7", "SMM8", "SMM9", "SMM10", "SMM11", 
                    "MM1", "MM2", "MM3", "MM4", "MM5", "MM6", "MM7"))

# Move the all.integrated object in the all_list to allow reuse of code
all_list <- list("all.integrated" = all.integrated)

# Switch default assay to "RNA" for plotting of DotPlots
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results. #NB: This will remove genes from the object
for (i in names (all_list)) {
  counts <- GetAssayData(all_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list[[i]] <- subset(all_list[[i]], features = rownames(counts))
}


# 4. Fix DE result lists to correct numbers of clusters and annotate them --------------------------
# 4.A DE testing - FindAllMarkers

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create a list and move DE results
FindAllMarkers_All_list <- list("all.integrated" = FindAllMarkers_All_list_5$all.integrated)

# Rename the clusters in the list to friendly names
FindAllMarkers_All_list$all.integrated = FindAllMarkers_All_list$all.integrated %>% 
  mutate(cluster = factor(cluster, levels = c(0:4), labels = cluster_names_Controls_2))

# Save out the DE results
write.xlsx(FindAllMarkers_All_list,
           "Output/DE/Controls_2/Final/FindAllMarkers_All_clusters.xlsx")

# Clean up the environment
rm(FindAllMarkers_All_list_6, FindAllMarkers_All_list_5)

# 4.A.ii DE testing - FindAllMarkers_Disease - Assess global differences between disease state -----

# Create a list and move DE results
FindAllMarkers_Disease_list <- FindAllMarkers_Disease_list_5

# Save out the DE results
write.xlsx(FindAllMarkers_Disease_list,
           "Output/DE/Controls_2/Final/FindAllMarkers_Disease.xlsx")

# Output Volcano plots for Disease restricted differences
EnhancedVolcano(FindAllMarkers_Disease_list$all.integrated, 
                lab = FindAllMarkers_Disease_list$all.integrated$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "Age matched controls vs. Disease",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/Controls_2/Final/",
                "HD vs. Disease", ".png"), scale = 2:1)

# Clean up the environment
rm(FindAllMarkers_Disease_list_6, FindAllMarkers_Disease_list_5)

# 4.A.iii DE testing - FindAllMarkers_Patient - Assess global differences between Patients ---------

# Create a list and move DE results
FindAllMarkers_Patient_list <- list("all.integrated" = FindAllMarkers_Patient_list_5$all.integrated)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list, 
           "Output/DE/Controls_2/Final/FindAllMarkers_Patient.xlsx")

# Clean up the environment
rm(FindAllMarkers_Patient_list_6, FindAllMarkers_Patient_list_5)

# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers Itemise differences between Clusters -------------------------------------------

# Create lists and move DE results
FindMarkers_All_list <- FindMarkers_All_list_5
# Rename list items to friendly names
names(FindMarkers_All_list) = cluster_names_Controls_2

# Save out the DE results
write.xlsx(FindMarkers_All_list, 
           "Output/DE/Controls_2/Final/FindMarkers_All.xlsx")

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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Final/Cluster/",
                                   names(FindMarkers_All_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_All_list_6, FindMarkers_All_list_5)

# 4.B.ii FindMarkers - Itemise differences between Disease -----------------------------------------

# Create a list and move DE results
FindMarkers_Disease_list <- FindMarkers_Disease_list_5

# Rename the items in the list
names(FindMarkers_Disease_list) <- cluster_names_Controls_2

# Save out the DE results
write.xlsx(FindMarkers_Disease_list, 
           "Output/DE/Controls_2/Final/FindMarkers_Disease.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Disease_list)) {
  temp = 
    (EnhancedVolcano(FindMarkers_Disease_list[[`i`]], 
                     lab = FindMarkers_Disease_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Disease_list)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Final/Disease/",
                                   names(FindMarkers_Disease_list)[i], ".png"), scale = 2:1)
}

# Clean up the environment
rm(FindMarkers_Disease_list_6, FindMarkers_Disease_list_5)


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Disease ----------------------------------------------

# Create a list and move DE results
FindConservedMarkers_Disease_list <- FindConservedMarkers_Disease_list_5

# Rename the items in the list
names(FindConservedMarkers_Disease_list) <- cluster_names_Controls_2

# Save out the DE results
write.xlsx(FindConservedMarkers_Disease_list, 
           "Output/DE/Controls_2/Final/FindConservedMarkers_Disease.xlsx")

# Clean up the environment
rm(FindConservedMarkers_Disease_list_6, FindConservedMarkers_Disease_list_5)

# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAP of cluster distribution of the all.integrated object
DimPlot(all.integrated, label = TRUE, split.by = "Disease", cols = ccolours) + 
  labs(title = "all.integrated", color = "Cluster") +
  theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Dimplots/Controls_2/Final/Dimplot_",
                "all.integrated", ".png"))


# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Create stacked bar plots of cluster distribution on a per disease level
# Gather data to plot
temp_Disease <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$Disease))
# Convert cell numbers to proportions based on disease
temp_Disease = temp_Disease %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by disease status
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 4:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster",
                    labels = rev(cluster_names_Controls_2)) +
  labs(x = "Patient cohort", y = "Cells", 
       title = "Controls_2 - BM cluster distribution by disease status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Final/Barplot_",
                "Controls_2 BM_clusters_by_disease", ".png"))

# Plot as proportion of total, split by disease
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 4:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_2)) + 
  labs(x = "Patient cohort", y = "Proportion (%)", 
       title = "Controls_2 - BM cluster distribution by disease status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Final/Barplot_",
                "Controls_2 BM_clusters_by_disease_proportions", ".png"))

# Create stacked bar plots of cluster distribution on a per sample level
# Gather data to plot
temp_all.integrated <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$orig.ident))
# Convert cell numbers to proportions based on disease
temp_all.integrated = temp_all.integrated %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_2),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  geom_text(aes(label = paste(rev(Freq))), position = position_stack(vjust = 0.5)) +
  labs(x = "Sample", y = "Cells", 
       title = "Controls 2 - BM cluster distribution by sample") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Final/Barplot_",
                "Controls_2_BM_clusters_by_sample", ".png"))

# Plot as proportion of total, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Proportion), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_2),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Controls 2 - BM cluster distribution by sample") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Final/Barplot_",
                "Controls_2_BM_clusters_by_sample_proportions", ".png"))

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by Disease
top_10_markers_all.integrated <- 
  FindAllMarkers_All_list$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all.integrated <- unique(top_10_markers_all.integrated$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all.integrated, split.by = "Disease", cols = c("orange", "turquoise", "green", "purple")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all.integrated") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_2/Final/Dotplot_",
                "split_by_disease", ".png"), scale = 2:1)

# 4.D.iv Heatmap of genes that define clusters based on level of expression ------------------------

# Delete the barcode metadata - not sure why this is required but OK
#temp = all.integrated
#temp$barcode = NULL
#
# Plot MulitHeatMap
#png("Output/Figures/Seurat/Controls_2/Final/Heatmap/Heatmap_all.integrated.png")
#plot_heatmap(temp,
#             n = 10,
#             markers = top_10_unique_all.integrated,
#             sort_var = c("Cluster","Disease"),
#             anno_var = c("Cluster", "Disease"),
#             anno_colors = list(ccolours,
#                                c("orange", "turquoise", "green", "purple")),
#             hm_limit = c(-2,0,2),
#             hm_colors = c("purple","black","yellow"))
#dev.off()

# 5. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Save the annotated Seurat objects for further analysis
saveRDS(all.integrated, "Data/R_out/Controls_2/Controls_2_complete.rds")

# Clean the environment and save out the DE analysis
rm(list = ls()[!ls() %in% c("FindAllMarkers_All_list", "FindAllMarkers_Disease_list",
                            "FindAllMarkers_Patient_list", "FindConservedMarkers_Disease_list",
                            "FindMarkers_All_list", "FindMarkers_Disease_list")])

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_2/DE_Analysis_Controls_2.rds")

print("#># Finished running '6. Controls dataset #2 - DE Analysis' script")