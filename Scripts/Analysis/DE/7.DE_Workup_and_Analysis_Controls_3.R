# Script information -------------------------------------------------------------------------------

# Title: Controls dataset #3 - DE Workup and analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will perform both the DE workup and analysis on the Controls 
# dataset #3 from Mogilenko et al., 2021: https://doi.org/10.1016/j.immuni.2020.11.005. 
# As the authors have provided detailed analysis we will attempt to recreate their findings. We will
# standardise the colours for clusters and ensure that plots output in a sensible manner.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '7. Controls dataset #3 - DE Workup and analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, openxlsx, multtest, metap, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed Control #1 dataset
integrated_list <- readRDS("Data/R_out/Controls_3/integrated_list.rds")

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")
annotations <- read.csv("Data/Public_data/annotations.csv") 

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

# Define cluster names for naming of lists - 5 clusters
cluster_names <- c("TN", "TTE", "TEM", "MAIT", "TCM")
cluster_names_Controls_3 <- c("TN", "TTE", "TEM", "MAIT", "TCM")

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")
# These need to be finalised
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
ccolours <- c("#00BA38", "#B79F00", "#F8766D", "#00BFC4", "#F564E3")


# 3. Pre-processing --------------------------------------------------------------------------------

# Copy the integrated_list to the all_list to allow copying of previous code and backup
all_list <- integrated_list

# Switch default assay to "RNA" for differential expression testing
DefaultAssay(all_list) <- "RNA"

# Remove non helpful genes from DE results. #NB: This will remove genes from the object
counts <- GetAssayData(all_list, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
all_list <- subset(all_list, features = rownames(counts))


# 4. Analysis of 5 cluster Seurat objects ----------------------------------------------------------
# 4.A DE testing - FindAllMarkers

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Run FindAllMarkers
FindAllMarkers_All_list <- FindAllMarkers(all_list)

# Annotate the results 
FindAllMarkers_All_list <-inner_join(x = FindAllMarkers_All_list, 
                                       y = annotations[, c("gene_name", "description")],
                                       by = c("gene" = "gene_name")) %>% unique()

# Rename the clusters in the list to friendly names
FindAllMarkers_All_list = FindAllMarkers_All_list %>% 
  mutate(cluster = factor(cluster, levels = c(0:4), labels = cluster_names))

# Save out the DE results
write.xlsx(FindAllMarkers_All_list,
           "Output/DE/Controls_3/Final/FindAllMarkers_All.xlsx")

# 4.A.ii DE testing - FindAllMarkers_Age - Assess global differences between disease state ---------
# Copy the all integrated objects into a new list and switch the identity class to "age"
all_Age <- all_list
Idents(all_Age) = "age"

# Run FindAllMarkers
FindAllMarkers_Age_list <- FindAllMarkers(all_Age)

# Annotate the results 
FindAllMarkers_Age_list <-inner_join(x = FindAllMarkers_Age_list, 
                                       y = annotations[, c("gene_name", "description")],
                                       by = c("gene" = "gene_name")) %>% unique()

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = old, pct.2 = young
FindAllMarkers_Age_list <- 
  FindAllMarkers_Age_list[!(FindAllMarkers_Age_list$cluster == "young"),]

# Save out the DE results
write.xlsx(FindAllMarkers_Age_list,
           "Output/DE/Controls_3/Final/FindAllMarkers_All.xlsx")


# 4.A.iii DE testing - FindAllMarkers_Patient - Assess global differences between Patients ---------

# Copy the all_list into a new list and switch the identity class to "orig.ident"
all_Patient <- all_list
Idents(all_Patient) = "orig.ident"

# Run FindAllMarkers
FindAllMarkers_Patient_list <- FindAllMarkers(all_Patient)

# Annotate the results 
FindAllMarkers_Patient_list <-inner_join(x = FindAllMarkers_Patient_list, 
                                           y = annotations[, c("gene_name", "description")],
                                           by = c("gene" = "gene_name")) %>% unique()

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list, 
           "Output/DE/Controls_3/Final/FindAllMarkers_Patient.xlsx")


# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 5 clusters
for (i in 1:5){
  temp <- FindMarkers(all_list, ident.1 = i-1)
  FindMarkers_All_list[[length(FindMarkers_All_list) +1]] <- 
    temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindMarkers_All_list) <- cluster_names

# Save out the DE results
write.xlsx(FindMarkers_All_list, 
           "Output/DE/Controls_3/Final/FindMarkers_All.xlsx")

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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_3/Final/Cluster/",
                                    names(FindMarkers_All_list)[i], ".png"), scale = 2:1)
}

# 4.B.ii FindMarkers - Itemise differences between Age ---------------------------------------------

# Create an empty list to store results in
FindMarkers_Age_list <- list()

# Run for loop to iterate across all clusters
for (i in 1:5){
  temp <- FindMarkers(all_list, 
                      ident.1 = "old", group.by = "age", subset.ident = i-1)
  FindMarkers_Age_list[[length(FindMarkers_Age_list) +1]] <- 
    temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindMarkers_Age_list) <- cluster_names

# Save out the DE results
write.xlsx(FindMarkers_Age_list, 
           "Output/DE/Controls_3/Final/FindMarkers_Age.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Age_list)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Age_list[[`i`]], 
                     lab = FindMarkers_Age_list[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Age_list)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_3/Final/Age/",
                                   names(FindMarkers_Age_list)[i], ".png"), scale = 2:1)
}


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Age -------------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Age_list = list()

# Run for loop to iterate across all clusters
for (i in 1:5) {
  temp <- FindConservedMarkers(all_list, 
                               ident.1 = i-1, grouping.var = "age")
  FindConservedMarkers_Age_list[[length(FindConservedMarkers_Age_list) +1]] <- 
    temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindConservedMarkers_Age_list) <- cluster_names

# Save out the DE results
write.xlsx(FindConservedMarkers_Age_list, 
           "Output/DE/Controls_3/Final/FindConservedMarkers_Age.xlsx")


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# Rename the integrated_list with annotated cluster names
integrated_list$Cluster = integrated_list$seurat_clusters
Idents(integrated_list) = "Cluster"
names(cluster_names) <- levels(integrated_list)
integrated_list <- RenameIdents(integrated_list, cluster_names)

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
DimPlot(integrated_list, split.by = "age", cols = ccolours, label = TRUE) + 
  labs(title = "Controls_3  - split by age", color = "Cluster") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Dimplots/Controls_3/Final/Dimplot_",
                "Controls_3 PB", ".png"))

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Create barplot of distribution by age
# Gather data to plot
temp_all.age <- as.data.frame(table(all_list$seurat_clusters, all_list$age))
# Convert cell numbers to proportions
temp_all.age <- temp_all.age %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells
ggplot(temp_all.age, 
       aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_3),
                    (rev(levels(all_list@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Cells", title = "Controls_3 PB") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_3/Final/Barplot_",
                "Controls_3_PB_clusters_by_age", ".png"))

# Plot as a proportion
ggplot(temp_all.age, 
       aes(fill = Var1, y = rev(Proportion), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_3),
                    (rev(levels(all_list@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Proportion (%)", title = "Controls_3 PB") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_3/Final/Barplot_",
                "Controls_3_PB_clusters_by_age_proportions", ".png"))

# Calculate distribution of clusters on a per sample level
# Gather data to plot
temp_all.integrated <- as.data.frame(table(all_list$seurat_clusters, all_list$orig.ident))
# Convert cell numbers to proportions
temp_all.integrated <- 
  temp_all.integrated %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster", 
                    labels = rev(cluster_names_Controls_3),
                    (rev(levels(all_list@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Cells", title = "Controls_3 PB") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_3/Final/Barplot_",
                "Controls_3_PB_clusters_by_sample", ".png"))

# Plot as proportion of total
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Proportion), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours), name = "Cluster",
                    labels = rev(cluster_names),
                    (rev(levels(all_list@meta.data$seurat_clusters)))) + 
  labs(x = "Patient cohort", y = "Proportion (%)", title = "Controls_3 PB") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_3/Final/Barplot_",
                "Controls_3 PB_clusters_by_sample_proportions", ".png"))


# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by Age
top_10_markers_all <- FindAllMarkers_All_list %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list, 
        features = top_10_unique_all, split.by = "age", 
        cols = c("green", "purple")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_3/Final/Dotplot_",
                "Age", ".png"), scale = 2:1)

# 4.D.iv Heatmap of genes that define clusters based on level of expression ------------------------

# Plot MulitHeatMap
#png("Output/Figures/Seurat/Controls_3/Heatmap/Heatmap_all.integrated.png")
#plot_heatmap(integrated_list,
#             n = 10,
#             markers = top_10_unique_all,
#             sort_var = c("Cluster","age"),
#             anno_var = c("Cluster", "age"),
#             anno_colors = list(ccolours,
#                                c("green", "purple")),
#             hm_limit = c(-2,0,2),
#             hm_colors = c("purple","black","yellow"))
#dev.off()

# 5. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Save the annotated control data
saveRDS(integrated_list, "Data/R_out/Controls_3/Controls_3_complete.rds")

# Clean up the environment
rm(list = 
     ls()[!ls() %in% c("annotations", "annotate_markers", "black_list_TCR", 
                       "FindAllMarkers_All_list", "FindAllMarkers_Age_list", 
                       "FindAllMarkers_Patient_list", "FindConservedMarkers_Age_list", 
                       "FindMarkers_All_list", "FindMarkers_Age_list")])

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_3/DE_Workup_and_Analysis_Controls_3.rds")

print("#># Finished running '7. Controls dataset #3 - DE Workup and analysis' script")